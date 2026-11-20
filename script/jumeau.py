"""MFM et tMFM_G. Adaptation Python des équations de Pathak et Battaglia."""
import numpy as np
from numba import njit
from scipy.signal import butter, sosfilt, sosfiltfilt

from commun import MODEL_COLORS, aligned_sc, figure, read_json, required, source, write_json, write_table
from dynamique import TARGETS, feature_set, mc_matrix, validation


@njit
def integrate(s, g, noise, noise_g, c, dt, fixed, param, ceiling, per_step):
    tau, gamma, jn, i0, w, a, b, d = fixed
    base, sigma, theta, volatility = param
    root = np.sqrt(dt)
    scale = sigma if per_step else sigma * root
    rates = np.empty_like(noise)
    for t in range(len(noise)):
        u = a * (jn*w*s + i0 + jn*g*(c @ s)) - b
        for i in range(len(s)):
            if abs(u[i]) < 1e-8:
                r = 1/d + u[i]/2
            elif -d*u[i] > 700:
                r = 0.0
            else:
                r = u[i] / (-np.expm1(-d*u[i]))
            s[i] += dt*(-s[i]/tau + (1-s[i])*gamma*r) + scale*noise[t,i]
            if not np.isfinite(s[i]) or abs(s[i]) > 1e5:
                raise ValueError("Intégration instable : vérifier dt, les unités et les bornes")
            rates[t,i] = r
        g = min(ceiling, max(base, g + theta*(base-g)*dt + volatility*root*noise_g[t]))
    return rates, g


def setup(cfg, model, names):
    import math
    required(model, ["sc_model","sc_model_labels","sc_model_weights","sc_normalization","neural_noise_convention"])
    files = {k: source(dict(model, dataset_root=cfg["dataset_root"]), k) for k in ("sc_model","sc_model_labels")}
    sc = aligned_sc(files["sc_model"], files["sc_model_labels"], names)
    if np.any(sc < 0) or not sc.max():
        raise ValueError("SC de simulation requise : poids non négatifs, matrice non nulle")
    if model["sc_normalization"] == "column_max":
        colmax = sc.max(axis=0)
        sc = sc / np.where(colmax>0, colmax, 1)
    elif model["sc_normalization"] == "global_max":
        sc = sc / sc.max()
    elif model["sc_normalization"] != "as_supplied":
        raise ValueError("Normalisation SC : column_max, global_max ou as_supplied")
    if model["neural_noise_convention"] not in ("diffusion_sqrt_dt","per_step"):
        raise ValueError("Convention de sigma : diffusion_sqrt_dt ou per_step")
    if type(model["global_signal_regression"]) is not bool or model["observation_band_hz"] is None:
        raise ValueError("Préciser global_signal_regression et observation_band_hz ([] pour aucun filtre supplémentaire)")
    tr, dt = cfg["tr_seconds"], model["dt_seconds"]
    if not 0 < dt < tr or not np.isclose(round(tr/dt)*dt, tr, rtol=0, atol=1e-8):
        raise ValueError("dt doit diviser le TR")
    if not 0 < model["lowpass_hz"] < 0.5/tr:
        raise ValueError("Le passe-bas de simulation doit être sous Nyquist au TR cible")
    if model["burn_seconds"] < 0 or not math.isfinite(model["burn_seconds"]):
        raise ValueError("Durée de transitoire non valide")
    return np.ascontiguousarray(sc), files


def simulate(sc, cfg, model, params, n, seed):
    tr, dt = cfg["tr_seconds"], model["dt_seconds"]
    stride = round(tr/dt)
    burn = round(model["burn_seconds"]/dt)
    total = burn + n*stride
    fixed_names = ("tau_s","gamma","J_N","I0","w","a","b","d")
    fixed = np.array([params.get(k,model["fixed"][k]) for k in fixed_names], dtype=float)
    par = np.array([params["G"],params["sigma"],params.get("theta",0),params.get("sigma_ou",0)], dtype=float)
    ceiling = float("inf") if model["ou_ceiling"] is None else float(model["ou_ceiling"])
    if not np.isfinite(fixed).all() or fixed[0]<=0 or fixed[-1]<=0 or not np.isfinite(par).all() or par.min()<0 or ceiling<par[0]:
        raise ValueError("Paramètres neuronaux/OU hors domaine")
    rng = np.random.default_rng(seed)
    s, g = np.full(len(sc),0.1), par[0]
    sos = butter(4, model["lowpass_hz"], fs=1/dt, output="sos")
    zi = np.zeros((len(sos), 2, len(sc)))
    chunks = []
    for start in range(0, total, 4096):
        count = min(4096,total-start)
        noise = rng.standard_normal((count,len(sc)))
        noise_g = rng.standard_normal(count)
        rates, g = integrate(s,g,noise,noise_g,sc,dt,fixed,par,ceiling,model["neural_noise_convention"]=="per_step")
        filtered, zi = sosfilt(sos, rates, axis=0, zi=zi)
        ix = np.arange(start,start+count) - burn
        chunks.append(filtered[(ix>=0) & (ix%stride==0)])
    y = np.concatenate(chunks)
    if model["global_signal_regression"]:
        centered = y - y.mean(axis=0)
        global_ts = centered.mean(axis=1)
        denom = global_ts @ global_ts
        if denom > 0:
            y = centered - np.outer(global_ts,(global_ts @ centered)/denom)
    band = model["observation_band_hz"]
    if band:
        lo, hi = band
        if not 0 <= lo < hi < 0.5/tr:
            raise ValueError("Bande d'observation incompatible avec le TR")
        filt = butter(2, [lo,hi] if lo else hi, btype="bandpass" if lo else "lowpass", fs=1/tr, output="sos")
        y = sosfiltfilt(filt, y, axis=0)
    return y


def fit(out, variant, cfg, model, sc, target, grid, n):
    keys = ["G","sigma"] + (["theta","sigma_ou"] if variant=="tmfm" else [])
    keys.extend(k for k in ("w","a","b","d") if model["bounds"].get(k) is not None)
    required(model["bounds"], keys)
    bounds = np.asarray([model["bounds"][k] for k in keys], dtype=float)
    if bounds.shape != (len(keys),2) or not np.isfinite(bounds).all() or not np.all(bounds[:,1] > bounds[:,0]) or bounds.min()<0:
        raise ValueError("Bornes attendues : [minimum, maximum] pour chaque paramètre")
    if model["ou_ceiling"] is not None and model["ou_ceiling"] < bounds[0,1]:
        raise ValueError("Le plafond de G doit couvrir sa plage ajustée")
    size, generations = model["population"], model["generations"]
    if size < 4 or generations < 1 or not model["optimizer_seeds"] or not model["fit_seeds"]:
        raise ValueError("Population, générations ou graines insuffisantes")
    if set(model["fit_seeds"]) & set(model["validation_seeds"]):
        raise ValueError("Séparer les graines d'ajustement et d'évaluation")
    trace, solutions = [], []
    for runseed in model["optimizer_seeds"]:
        rng = np.random.default_rng(runseed)
        population = rng.uniform(bounds[:,0], bounds[:,1], size=(size, len(keys)))
        cache = {}
        for generation in range(generations):
            scores = []
            for individual in population:
                key = tuple(individual.tolist())
                if key not in cache:
                    prediction = []
                    try:
                        for seed in model["fit_seeds"]:
                            y = simulate(sc,cfg,model,dict(zip(keys,key)),n,seed)
                            prediction.append(feature_set(y[grid],grid,cfg["tr_seconds"])[0])
                        cache[key] = np.mean([np.linalg.norm(v-target) for v in prediction])
                    except ValueError as exc:
                        if not str(exc).startswith(("Intégration instable","Corrélation indéfinie")):
                            raise
                        print(f"Candidat écarté, {dict(zip(keys,key))} : {exc}",flush=True)
                        cache[key] = np.inf
                scores.append(cache[key])
            if not np.isfinite(scores).any():
                raise ValueError("Population entière non estimable : revoir les bornes et le contrôle de simulation")
            order = np.argsort(scores)
            best = population[order[0]].copy()
            loss = float(scores[order[0]])
            trace.append([runseed,generation+1,loss,*best])
            write_table(out / "data/trace.csv", ["optimizer_seed","generation","error",*keys], trace)
            parents = population[order[:max(1,min(np.isfinite(scores).sum(),size//2))]]
            offspring = []
            for _ in range(size-1):
                p,q = parents[rng.integers(len(parents),size=2)]
                child = np.where(rng.random(len(keys))<0.5,p,q)
                mutate = rng.random(len(keys)) < 1/len(keys)
                child += mutate*rng.normal(0,0.1,len(keys))*(bounds[:,1]-bounds[:,0])
                offspring.append(np.clip(child,bounds[:,0],bounds[:,1]))
            population = np.vstack([best,offspring])
        solutions.append(dict(optimizer_seed=runseed,error=loss,parameters=dict(zip(keys,best.tolist()))))
    write_json(out / "data/solutions.json", solutions)
    write_json(out / "data/best.json", min(solutions,key=lambda s:s["error"]))
    def draw(fig, ax):
        for i, seed in enumerate(model["optimizer_seeds"]):
            v = [r for r in trace if r[0]==seed]
            ax.plot([r[1] for r in v],[r[2] for r in v],label=str(seed), color=MODEL_COLORS[variant],
                    linestyle=("-", "--", ":", "-.")[i % 4], marker=("", "o", "s", "^")[(i // 4) % 4])
        ax.set(xlabel="Génération",ylabel="Erreur euclidienne moyenne",title=variant+" · ajustements")
        ax.legend(title="Graine")
    figure(out, "ajustement", draw)


def evaluate(out, cfg, model, sc_model, sc_obs, fits, empirical, grid, n):
    seeds = model["validation_seeds"]
    if not seeds or set(seeds) & set(model["fit_seeds"]):
        raise ValueError("Renseigner des graines de validation indépendantes")
    rows = []
    for variant, path in fits.items():
        params = read_json(path)["parameters"]
        for seed in seeds:
            y = simulate(sc_model,cfg,model,params,n,seed)[grid]
            target, stream, *_ = feature_set(y,grid,cfg["tr_seconds"])
            vals, _ = validation(y, sc_obs, stream, mc_matrix(stream))
            vals.update(dict(zip(TARGETS, target)))
            for key, v in vals.items():
                rows.append([variant,seed,key,empirical[key],v,v-empirical[key]])
    write_table(out / "data/evaluation.csv", ["model","seed","feature","empirical","simulated","difference"],rows)
    def draw(fig, ax):
        keys = [k for k in empirical if k not in TARGETS]
        ax.plot(keys,[empirical[k] for k in keys],"ko-",label="Empirique")
        for modelname in fits:
            pred = [[r[4] for r in rows if r[0]==modelname and r[2]==k] for k in keys]
            ax.plot(keys,np.mean(pred,axis=1),".-",label=modelname, color=MODEL_COLORS[modelname])
        ax.tick_params(axis="x",rotation=45)
        ax.set(ylabel="Quantile",title="Mesures réservées à l’évaluation")
        ax.legend()
    figure(out,"validation",draw)

