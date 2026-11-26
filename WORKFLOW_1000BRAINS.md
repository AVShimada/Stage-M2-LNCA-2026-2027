# 1000BRAINS · du premier sujet V4 au jumeau mécanistique

## Relance complète depuis un dossier `output` neuf

La relance reconstruit les descriptifs démographiques et cognitifs des **409 participants**, puis les analyses empiriques, les graphes et le jumeau technique de **sub-0001 / ses-1 / V4**. Les comparaisons V1/V4 utilisent la même acquisition et le même atlas Schaefer 100. Le lanceur attribue un identifiant UTC commun et relit les résultats correspondants à chaque étape.

Après avoir renommé l’ancien dossier `output` à la racine du projet, exécuter dans une allocation de calcul sur le cluster :

```bash
cd /home/aurelio/Server/Aurelio_Stage_M2 &&
unzip -o maj_architecture_1000brains.zip &&
"$HOME/.venvs/lnca-1000brains/bin/python" -m pip install -r env/requirements-1000brains.txt &&
"$HOME/.venvs/lnca-1000brains/bin/python" script/relancer_workflow.py
```

| Dans `output/` | Ordre de calcul |
|---|---|
| `01_population/` | 01 métadonnées → 02 scores cognitifs de visite 1 → 03 audit des tests manquants → 04 deuxième scan déclaré → 05 tests attendus complets selon l’âge → 06 descriptifs des 409. |
| `02_pilote/` | 01 contrôle des entrées → 02 accord FC → 03 aperçu SC et FC V1/V4 → 04 SC et FC statiques → 05 dynamique M1 → 06 FCD V1/V4 → 07 cibles 2026 → 08 communautés et hubs → 09 FC et distance → 10 rich club → 11 modules de méta-connectivité. |
| `03_jumeau/` | 01 simulation technique du MFM → 02 comparaison du BOLD simulé aux données empiriques. Les figures SC, FC, FCD, vitesses et cibles sont rassemblées sous `02_comparaison/<RUN>_sub-0001_ses-1_V4/figures/`. |
| `04_exports/<RUN>/` | Configuration effective du pilote, journal `relance.log` et archive ZIP des figures, README et manifestes. |

Le premier bloc prépare le score sur sa population de référence puis le relie à la sélection finale. Le filtre applique la règle AKT à partir de 55 ans. Les descriptifs doivent retrouver **409 sujets et 409 scores** pour poursuivre. Dans le bloc imagerie, la fenêtre de visualisation de la FCD est de **10 volumes**, le pas d’un volume et l’échelle de 0 à 0,5. Les dix cibles numériques conservent le protocole 2026.

Le lanceur privilégie `config/pilote_ordre_M1_20261120T150810Z.json`. Si cette fiche est absente, il lit `config/pilote.json`. Il complète les chemins anatomiques vides à partir de `atlas/labels_M1.txt` et `atlas/coordonnees_M1.csv` dans la copie du lancement. L’ordre commun des régions hérité du M1 reste une hypothèse anatomique. La partition en modules porte sur une matrice de **4 950 × 4 950** liens. Prévoir la mémoire et le temps de calcul correspondants dans l’allocation.

Les ajustements MFM/tMFM, l’imagerie de cohorte et les deltas entre visites correspondent aux phases suivantes. Chaque dossier de ces phases est créé lors de son premier calcul. Les chemins sont définis dans `script/commun.py`. Le lanceur régénère `architecture_projet.txt` en fin de calcul.

## 1. Entrées et conventions du pilote

### Fiche à remplir dans `config/pilote.json`

Le pilote est prérempli pour **sub-0001, ses-1, Schaefer 100 / 7 réseaux** à partir de l’inventaire communiqué depuis le cluster. La série BOLD V4, sa fiche TSdim, la FC Fisher z V4 et la SC SIFT2 CM_log figurent dans la configuration. Le même enregistrement dispose également des fichiers fonctionnels V1. Le fichier de pourcentage de volumes portant les mentions 400P et V2 concerne une autre combinaison d’atlas et de version. Le TR, la chronologie des volumes censurés et l’ordre anatomique restent à documenter.

Avant de renseigner ces métadonnées, exécuter sur le cluster :

```bash
python3 script/pipeline.py entrees --config config/pilote.json
```

Cette lecture affiche les dimensions sur disque et le début de TSdim. Elle repère l’axe des 100 régions si les dimensions le permettent, contrôle la symétrie des matrices et compare la FC de Pearson des échantillons livrés à la tangente hyperbolique de la FC Fisher z fournie, sur les 4 950 liens. Un écart mérite d’examiner la sélection des volumes, l’ordre des régions et la méthode de calcul. Même un accord numérique ne renseigne ni le TR ni les intervalles entre volumes conservés.

**Contrôle du 19 novembre 2026 sur le cluster.** Pour `sub-0001/ses-1` en V4, la série livrée et sa fiche TSdim indiquent 296 échantillons × 100 régions. La SC et la FC sont de taille 100 × 100. L’orientation confirmée est `time_regions`. Sur les 4 950 liens, la corrélation entre la FC de Pearson recalculée et la FC fournie reconvertie en r vaut **0,999999999999929**. La FC recalculée est supérieure sur **4 809 liens (97,15 %)**. L’écart médian signé vaut **+0,00788728** et l’écart absolu maximal **0,0187704**. Ces chiffres proviennent des sorties transmises depuis le cluster. Ils indiquent un profil de connectivité pratiquement identique avec une différence systématique de valeurs. Son origine, le TR, la censure et l’ordre des labels SC restent à documenter.

Le lancement `20261119T141004Z_sub-0001_ses-1_V4` a mesuré la relation **FC recalculée = 1,020751114 × FC fournie − 5,5343 × 10⁻¹⁰**, sur l’échelle de Pearson. L’écart résiduel absolu maximal est de **3,84087 × 10⁻⁷**. Pour cette acquisition, le facteur multiplicatif décrit donc presque entièrement la différence entre les matrices. Ce résultat porte sur la FC statique. Dans l’ancienne arborescence, la figure produite sur le cluster était `output/00_accord_fc/20261119T141004Z_sub-0001_ses-1_V4/figures/accord_fc.png` ; le même dossier contient son PDF. L’environnement Python utilisé avec succès pour le tracé s’active avec `source ~/.venvs/lnca-1000brains/bin/activate`.

**Durées variables.** Le pilote V4 conserve tous les volumes utilisables. Le paramètre `m1.max_volumes_for_speed` vaut `null` et la vitesse multi-échelles utilise elle aussi toute la série disponible. Les fenêtres restent à l’intérieur des segments continus. La durée simulée du jumeau suit la grille temporelle du sujet et les mêmes volumes sont sélectionnés pour l’extraction des cibles.

**Aperçu SC et FC statique.** Lancer `bash script/apercu_sc_fc.sh` pour tracer la SC SIFT2 `CM_log` et les FC Fisher z V1 et V4 livrées pour `sub-0001 / ses-1`. Les cartes reprennent la figure M1 non triée (diapo 33) avec la palette `turbo`, des cases carrées, le nœud 0 en haut et la diagonale à zéro à l’affichage. Les bornes de référence proviennent des matrices moyennes archivées des deux groupes du M1. Elles sont d’environ −0,401 à 1,293 pour la FC et −0,083 à 4,608 pour la SC. La figure `fc_v4_adaptee` utilise les 2e et 98e percentiles des 4 950 liens V4 pour mieux distinguer les blocs. Les figures `fc_v1` et `fc_v4` à échelle commune servent à comparer le niveau de FC. Les triangles des barres de couleurs signalent les valeurs hors échelle et le manifeste en donne les effectifs. Les figures sont écrites dans `output/02_pilote/03_apercu_sc_fc/`. La comparaison anatomique des liens SC et FC à l’étape `statique` utilise les listes de régions dans l’ordre de chaque fichier.

**Première FCD du pilote.** Le TR de 2 s provient des scripts V1. L’utilisateur indique que les 296 lignes V4 sont normalement consécutives. Pour ce premier essai, `censor_mode = "none"` signifie que l’analyse conserve toutes les lignes livrées, dans leur ordre et sans exclusion supplémentaire. Le protocole de débruitage reste à documenter. Cette hypothèse est propre au pilote. Une fenêtre de 60 s contient 30 volumes ; le pas de 2 s donne 267 fenêtres et une FCD de taille 267 × 267. Les 10 quantiles sont calculables sur le BOLD seul. L’alignement des labels SC et BOLD sera nécessaire pour le couplage SC–FC et le modèle structurel.

**Comparaison à 10 volumes.** Ce bloc reprend la série V4 complète de 296 volumes de `sub-0001`, avec des fenêtres de 10 volumes et un pas d’un volume. La carte utilise la palette `jet` de l’article et les bornes du M1, de 0 à 0,5, pour rendre visibles les différences hors diagonale. Le triangle supérieur de la barre de couleurs signale les valeurs au-dessus de 0,5. On attend une matrice de 287 × 287 fenêtres.

L’étape `06_fcd_v1_v4` enregistre les deux matrices et leurs figures dans le même lancement. Pour relancer uniquement cette comparaison :

```bash
bash script/comparer_fcd_v1_v4.sh config/pilote.json
```

### Couleurs des figures

| Figure | Convention d’affichage |
|---|---|
| FCD du M1 | Diapos 46 à 48. Palette `turbo` sur fond noir, bornes 0 et 0,5. |
| FCD du jumeau, empirique et simulée | Fenêtres de 10 volumes avancées d’un volume, palette `jet`, échelle de 0 à 0,5, fenêtre 0 en haut à gauche. |
| FC du jumeau, empirique et simulée | Référence `fc_v4_adaptee`. Les deux cartes utilisent `turbo` et la même échelle en z de Fisher, du 2e au 98e percentile des 4 950 liens empiriques. La carte de différence utilise une échelle divergente centrée sur zéro. |
| FC et SC du M1 | Palette `turbo`. La FC statique individuelle utilise les z de Fisher, le nœud 0 en haut et les percentiles 2 et 98 des liens comme bornes. Les moyennes des groupes partagent les bornes par type de connexion. Les différences entre versions ont des bornes symétriques autour de zéro. |
| Méta-connectivité | Palette `turbo`. Pour les moyennes par tertile et les modules, bornes 0 et 0,1 (diapos 49 et 50). Pour l’exemple individuel, le maximum est le 99e percentile des liens affichés. |
| Vitesses du M1 | Densités par noyau gaussien. Empirique en bleu foncé continu, phases randomisées en vert tireté, fenêtres mélangées en rouge pointillé (diapo 51). |
| Comparaison des jumeaux | Empirique en noir, MFM en bleu et tMFM à couplage variable en rouge orangé (figure 4B). |
| Âge et graphes | Âge continu en `viridis` (diapos 53 et 55). Force en cyan et participation en orange (diapos 42 et 43). Les courbes individuelles de distance et de rich club utilisent un trait blanc sur fond noir, le code de couleur par groupe d’âge étant réservé aux groupes. |

Les triangles aux extrémités d’une barre de couleurs signalent les valeurs qui dépassent les bornes affichées.

L’étape `02_pilote/07_cibles_2026` affiche la FCD de 10 volumes dans `figures/fcd.png` et conserve la matrice affichée dans `data/fcd_10_volumes.npy`. Le calcul des dix cibles utilise la FCD à 60 s enregistrée dans `data/fcd.npy` pour χ et les vitesses à 55–65 s pour ν. Les figures de densité lissent les distributions. Pour retracer les cartes de SC/FC, de dynamique M1 et de cibles d’un pilote enregistré dans cette arborescence, reprendre son identifiant de lancement dans la commande suivante.

```bash
cd /home/aurelio/Server/Aurelio_Stage_M2
"$HOME/.venvs/lnca-1000brains/bin/python" script/rendu_pilote.py --run "$RUN"
```

Le seuil de 280 volumes est rapporté de mémoire pour une analyse M1. Plusieurs scripts M1 archivés indiquent 290. Avant la cohorte, inventorier par session le nombre de volumes livrés, le nombre retenu après censure, le TR et le nombre de fenêtres exploitables. Si 280 devient un seuil d’inclusion, documenter à quel stade il s’applique. Un sujet au-dessus du seuil peut conserver toute sa série. Pour tester l’effet de la durée, comparer sur les mêmes sujets les résultats complets et les résultats à durée commune (p. ex. 280 volumes contigus si les TR sont identiques). Pour les vitesses M1, fixer `max_volumes_for_speed` à 280 limite le calcul aux indices de la grille initiale allant de 0 à 279. Avec des volumes censurés, cet intervalle peut contenir moins de 280 observations. Une comparaison à durée commune de la FC, de la FCD et du jumeau nécessite une sélection temporelle cohérente pour toutes ces analyses.

Pour mesurer séparément la pente et l’ordonnée à l’origine, lancer le contrôle suivant sur le cluster. Il ajuste `r_recalculé = pente × r_fourni + intercept` et affiche le résidu absolu maximal. La figure représente l’écart signé en fonction de la FC fournie. Une droite horizontale correspond à un décalage constant, une droite inclinée à une différence qui varie avec le poids du lien. Les sorties PNG, PDF, CSV et JSON sont enregistrées sous `output/02_pilote/02_accord_fc/`. Le calcul utilise les échantillons livrés du pilote, avant toute sélection supplémentaire de volumes.

```bash
python3 script/pipeline.py accord-fc --config config/pilote.json --run "$(date -u +%Y%m%dT%H%M%SZ)"
```

Le calcul numérique réalisé précédemment sur les mêmes 4 950 liens hors diagonale peut être reproduit depuis la racine du projet.

```bash
PYTHONPATH=script python3 - <<'PY'
import numpy as np
from commun import EDGE, array_file, read_json, source
from dynamique import fc

c = read_json("config/pilote.json")
x = array_file(source(c, "bold"))  # pilote confirmé : temps × régions
r_calc = fc(x)[EDGE]
r_fourni = np.tanh(array_file(source(c, "fc_reference_z"))[EDGE])
delta = r_calc - r_fourni
print("Corrélation entre les deux profils de liens :", np.corrcoef(r_calc, r_fourni)[0, 1])
print("Écart médian signé (recalculée − fournie) :", np.median(delta))
print("Liens où la FC recalculée est supérieure (%) :", 100 * np.mean(delta > 0))
PY
```

| Champ | Valeur attendue |
|---|---|
| `subject`, `session` | Identifiants exacts du même individu et de la même acquisition que les fichiers. |
| `version`, `atlas` | Déjà renseignés : V4 et Schaefer 100. |
| `dataset_root`, `bold` | Racine absolue du dataset et chemin exact de la **série temporelle BOLD régionale** V4. Un fichier de FC déjà calculé ne permet pas de reconstruire la FCD. |
| `sc`, `sc_weights` | Matrice structurelle correspondante et convention des poids. Le M1 utilisait notamment SIFT2 CM_log. Consigner celle réellement livrée. |
| `bold_orientation`, `tr_seconds` | `time_regions` ou `regions_time`, puis TR en secondes d’après les métadonnées. |
| `labels_bold`, `labels_sc` | Deux listes dans l’ordre propre à chaque fichier, 100 labels distincts, un par ligne. Le lecteur réaligne la SC sur l’ordre BOLD. Si la documentation garantit le même ordre que le M1, les deux chemins peuvent pointer vers `atlas/labels_M1.txt` avec son chemin absolu sur le cluster. |
| `censor_mode` | `none` pour une grille complète et continue sans volume à exclure, `keep_mask` si une ligne BOLD existe pour chaque volume et qu’un masque indique lesquels garder, `retained_indices` si le fichier BOLD a déjà été raccourci. |
| `censor_file`, `original_volume_count` | Pour un masque : une colonne de 0 et 1, 1 = garder. Pour une série raccourcie : un indice d’origine **en base zéro** par ligne conservée, plus le nombre total de volumes avant suppression. Laisser ces champs vides lorsqu’ils ne s’appliquent pas. |
| `preprocessing_reference` | Référence du README ou du rapport de débruitage dès qu’il est disponible. Relever mouvement, régression du signal global, filtrage, interpolation et censure. |
| `coordinates_csv` | Seulement pour FC–distance. Chemin du CSV avec `roi` ou `label`, et `x,y,z` en mm. `atlas/coordonnees_M1.csv` est fourni comme référence, avec association par nom de région. |

Les chemins d’entrée relatifs sont résolus depuis `dataset_root` ; un chemin absolu peut aussi être donné. Les lecteurs acceptent des matrices numériques sans en-tête en TXT/TSV (espaces ou tabulations), CSV (virgules) et NPY. Si la livraison utilise MAT, NIfTI, des en-têtes ou des index différents, adapter d’abord le lecteur à partir d’un exemple. Conserver la chronologie réelle des volumes lors de la conversion. Les visites déclarées proviennent de `Demographic_data.csv`.

```bash
python3 script/pipeline.py champs --config config/pilote.json
```

Cette commande liste les champs encore vides, y compris les champs conditionnels. Les paramètres inconnus sont représentés par `""` ou `null`. Aucun TR, sujet ou fichier V4 fictif n’est proposé.

## 2. Lancer le pilote V4, puis regarder ses figures

Sur un nœud de calcul, à la racine du projet :

```bash
CFG=config/pilote_ordre_M1_20261120T150810Z.json
RUN=$(date -u +%Y%m%dT%H%M%SZ)
echo "$RUN"
bash script/lancer_1000brains.sh pilote "$CFG" "$RUN"
bash script/lancer_1000brains.sh figures "$CFG" "$RUN"
```

Conserver la valeur affichée de `RUN` dans le cahier de laboratoire pour les jobs suivants. Le raccourci `pilote` exécute les quatre premières lignes du tableau **dans cet ordre**. Une étape peut aussi être appelée seule avec `bash script/lancer_1000brains.sh ETAPE "$CFG" "$RUN"`.

| Ordre et commande | Données produites | Contrôle à faire |
|---|---|---|
| 01 · `controle` | `quality.json`, ordre des régions, indices et temps retenus, figure de trois séries | Vérifier le nombre initial de volumes, le nombre retenu, les segments continus, le TR et la convention de SC. |
| 04 · `statique` | SC alignée, FC Pearson, FC en z de Fisher, table des 4 950 liens, métriques globales et figures | Examiner la symétrie, la distribution des poids et les correspondances anatomiques. Utiliser la moyenne **FC z** pour comparer aux archives M1. |
| 05 · `m1` | FC successives, FCD et MC à 20 s, vitesses à fenêtre fixe et multi-échelles, substituts, Vtyp et largeurs omises | Inspecter FCD, distribution des vitesses et nombre de paires par largeur. Les périodes censurées interrompent les fenêtres. |
| 07 · `cibles` | FCD à 60 s, vitesses de 55–65 s, dix cibles dans `targets.json`, FC/λ/MC de validation | Lire le nombre de valeurs de χ et de ν, les secondes effectivement réalisées et les distributions avant d’ajuster un modèle. |

**Définitions qui changent l’interprétation.** Pour deux fenêtres, χ est la corrélation entre leurs vecteurs de 4 950 liens FC. La distribution de χ prend le triangle supérieur strict de la FCD, y compris les paires de fenêtres qui se recouvrent. Pour ν, on compare deux fenêtres successives de même longueur qui se touchent sans recouvrement : `ν = 1 − corr(FC₁, FC₂)`. Les valeurs peuvent dépasser 1. Les cinq percentiles 5, 25, 50, 75 et 95 de chaque distribution donnent les dix cibles. **σ est estimé ensuite pendant l’ajustement**, en même temps que G.

- **Profil M1.** FCD/MC sur 20 s, pas d’un TR. Pour Vtyp moyen (15–60 s) et long (60–210 s), le réglage livré `pooled_sliding` suit le calcul de `dFC_speed.py` sur toute la durée utilisable : pas d’un TR, comparaison à une largeur de fenêtre, tailles espacées de 10 TR, médiane de toutes les vitesses réunies. `median_nonoverlap` reproduit l’autre convention de `dFC_speed4.py` : fenêtres disjointes et médiane des médianes par largeur. Consulter `omissions.json` si le scan est court. La provenance exacte de la variante des CSV M1 reste à confirmer avant une comparaison numérique stricte.
- **Substituts.** Même rotation de phase par fréquence pour toutes les régions, appliquée dans chaque segment continu. Le second substitut permute l’ordre des FC fenêtrées. Une réalisation de chaque type est prévue pour le premier contrôle visuel. Le nombre est réglable et enregistré.
- **Profil 2026.** Fenêtres de 60 s avancées de 2 s, arrondies au TR réalisable. Pour ν, toutes les largeurs entières en TR comprises entre 55 et 65 s sont regroupées. Cette grille est notre adaptation au TR livré. L’origine des fenêtres repart au début de chaque segment continu. `windows.json` et `speeds.csv` indiquent les pas, largeurs et effectifs réels. Les mesures de validation portent sur les entrées de FC (ρ), le couplage SC–FC de chaque fenêtre (λ) et les paires de liens de la MC. La représentation compacte de MC utilise 4 950 liens uniques pour éviter leur duplication.

### Figures complémentaires du pilote

| Ordre | Commande courte | Prérequis et sortie |
|---|---|---|
| 08 | `graphes` | Après `statique`. Louvain répété, consensus FC/SC, force et participation des régions, efficacité et clustering. Les arêtes positives sont utilisées pour les graphes, les poids descriptifs sont conservés à l’étape 04. |
| 09 | `distance` | Après `statique`, avec les coordonnées renseignées. FC positive par classe de distance en mm. |
| 10 | `rich-club` | Après `statique`. Courbe FC pondérée, 100 graphes nuls par défaut, graines et valeurs brutes. Un seuil sans coefficient défini donne une case vide ou NaN. |
| 11 | `modules-mc` | Après `m1`. Louvain sur la MC positive, module de chaque lien, matrice triée. La MC complète a 24 502 500 entrées et l’algorithme en conserve plusieurs copies. |

Ces quatre commandes peuvent être lancées indépendamment dès que leurs entrées existent. Pour la comparaison avec V1, dupliquer `config/pilote.json` en `config/pilote_v1.json`, passer les chemins `bold`, `ts_dimensions` et `fc_reference_z` sur les trois fichiers V1 du même répertoire, puis renseigner les métadonnées, les éventuels masques et `version = "V1"`. Puis :

```bash
RUN_V1=$(date -u +%Y%m%dT%H%M%SZ)
bash script/lancer_1000brains.sh pilote config/pilote_v1.json "$RUN_V1"
bash script/lancer_1000brains.sh comparer "$CFG" "$RUN" --reference config/pilote_v1.json --reference-run "$RUN_V1"
```

Le résultat est **V4 moins V1, sur la même acquisition**. Si les volumes retenus diffèrent, préparer deux configurations dédiées sur leur intersection et relancer les deux pilotes pour cette comparaison. Documenter d’abord les effectifs exclus par chaque version. Une différence de Vtyp chez un individu ne suffit pas à attribuer un effet au vieillissement ou au débruitage.

## 3. Construire puis évaluer le jumeau de ce sujet

**Retracer une simulation existante du pilote V4.** La commande suivante retrouve le dernier contrôle terminé de `sub-0001 / ses-1 / V4` sous `output/03_jumeau/`. Pour choisir une simulation précise, ajouter `--source-run` suivi du chemin de son dossier. Les figures sont écrites sous `output/03_jumeau/02_comparaison/`. La relance complète fournit ce chemin source explicitement.

```bash
cd /home/aurelio/Server/Aurelio_Stage_M2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$HOME/.venvs/lnca-1000brains/bin/python" script/audits/comparer_jumeau.py
```

- `figures/02_fc_statique.png` compare la FC Fisher fournie pour le pilote V4, la FC simulée et leur différence en z. Les bornes des deux premières cartes viennent du fichier V4 utilisé pour `fc_v4_adaptee`. Lorsqu’un masque supplémentaire est appliqué au BOLD, la référence est recalculée sur les seuls volumes retenus. Dans `data/comparaison.json`, la corrélation de profil et les erreurs Pearson comparent les deux séries BOLD sur cette même grille.
- `figures/01_sc_du_sujet.png` affiche la SC individuelle SIFT2 `CM_log` livrée, sur l’échelle du M1, puis les poids injectés dans le modèle. Pour le contrôle par défaut, ces derniers sont divisés par le maximum de chaque colonne. Si le modèle utilise les poids bruts, la légende indique la normalisation globale appliquée uniquement à l’affichage.
- `figures/03_fcd_10_volumes.png` compare les deux FCD à 10 volumes. Pour 296 volumes consécutifs, chaque carte fait 287 × 287. Les matrices visualisées sont `data/fcd_empirique_10vol.npy` et `data/fcd_simulee_10vol.npy`. Les dix quantiles d’ajustement conservent leur profil 2026.

Compléter `config/jumeau.json` après inspection du pilote. Le fichier individuel SIFT2 `dseg_CM.txt` sans suffixe log est proposé pour `sc_model`. Vérifier ses unités, l’ordre des labels et la normalisation retenue (`global_max`, `column_max` ou `as_supplied`). Les échelles de G dépendent de cette décision.

| Réglage à fixer | Utilisation |
|---|---|
| `control_parameters.G`, `control_parameters.sigma` | Un point de contrôle pour vérifier la stabilité, le temps de calcul et la production du BOLD simulé. |
| `bounds.G`, `bounds.sigma` | Bornes numériques [min, max] à établir par balayage sur la SC Schaefer à l’échelle choisie. |
| `bounds.theta`, `bounds.sigma_ou` | Bornes du taux de rappel (s⁻¹) et de la volatilité du processus OU de G. Sa moyenne de rappel est le G de base. Le plancher est G ; `ou_ceiling: null` laisse le plafond ouvert selon les équations de l’article. |
| `neural_noise_convention` | Choisir `diffusion_sqrt_dt` pour une amplitude de diffusion continue ou `per_step` pour une amplitude par pas comme dans le balayage archivé. Une même valeur numérique de σ change de sens entre ces conventions. |
| `global_signal_regression`, `observation_band_hz` | Renseigner `true` ou `false` et la bande [bas, haut] en Hz, ou `[]` si aucun filtrage supplémentaire n’est indiqué. Aligner l’observation sur la documentation empirique. Une méthode de débruitage V4 plus complexe demandera une adaptation spécifique. |

```bash
bash script/lancer_1000brains.sh simuler "$CFG" "$RUN" --model config/jumeau.json
# Inspecter le BOLD simulé, la FCD et le temps de calcul avant les ajustements.
bash script/lancer_1000brains.sh mfm     "$CFG" "$RUN" --model config/jumeau.json
bash script/lancer_1000brains.sh tmfm    "$CFG" "$RUN" --model config/jumeau.json
bash script/lancer_1000brains.sh evaluer "$CFG" "$RUN" --model config/jumeau.json
bash script/lancer_1000brains.sh figures "$CFG" "$RUN"
```

Les deux ajustements dépendent du contrôle de simulation et peuvent occuper des jobs distincts. Le MFM estime G et σ. Le tMFM_G ajoute θ et σ_OU ; w, a, b et d sont fixés aux valeurs du modèle de référence livré, sauf ajout explicite de bornes. Chaque candidat est simulé sur la durée du scan d’origine, après 50 s de transitoire. Les volumes retenus et la fonction d’extraction des cibles sont identiques pour données et simulations. Le critère est la distance euclidienne entre les dix quantiles, moyennée sur les graines d’ajustement.

**Sorties à lire.** `trace.csv` montre la convergence par génération, `solutions.json` la dispersion entre ajustements et `best.json` les paramètres de la meilleure solution. `evaluation.csv` confronte les modèles sur des graines distinctes, y compris sur ρ, λ et MC réservés à l’évaluation. Un paramètre qui atteint une borne ou varie fortement entre répétitions demande une vérification du balayage avant interprétation. Cette évaluation utilise toujours la même acquisition empirique.

**Lien avec l’article.** L’archive versionnée a été inspectée (MD5 `c6d033a94c66bf50d2c53d0d4fe7262c`). Elle fournit `WongWangFluc.m`, un balayage et dfcWalk. Le script d’ajustement individuel complet et ses bornes ne figurent pas dans cette version. Le balayage exécute un pas de 0,0001 s, applique le bruit neuronal par pas et impose un plafond à G de 0,5. Le texte décrit un pas de 0,001 s et une troncature inférieure. Notre implémentation utilise par défaut le pas du texte, un passe-bas Butterworth à 0,2 Hz avant échantillonnage, la durée individuelle et un algorithme génétique Python avec élitisme, croisement et mutation. Le code archivé utilise un filtre FFT. La population de 50 et le maximum de 100 générations viennent de l’article. Le réglage complet est coûteux et l’implémentation fournie évalue les candidats séquentiellement. Mesurer le coût du contrôle avant de dimensionner le job. Une réplication numérique exacte du pipeline MATLAB reste à vérifier. [Article](https://doi.org/10.1371/journal.pbio.3003916) · [Code versionné](https://doi.org/10.5281/zenodo.20798973)

## 4. Métadonnées de population et future imagerie de cohorte

### Préparer les métadonnées sources

Les trois fichiers récupérés sont sous `1kBrains/metadonnees/`. `metadonnees.py` prépare la démographie et un inventaire des visites à partir de `Demographic_data.csv`. Il relève aussi les noms, dimensions et types des variables de `cog_scores.mat`. Le score M1 est calculé directement par `cognition.py` à partir de la feuille `Sheet1` du classeur `1000Brains_NP_VarOI_V1_Pseudo.xlsx`. Il utilise la visite 1, les 19 tests orientés, la règle AKT à partir de 55 ans et la moyenne des scores standardisés disponibles. Comme au M1, la standardisation utilise les lignes dont l’âge est renseigné. Le suffixe du classeur désigne la version des données cognitives.

L’inventaire MATLAB permettra de vérifier la définition des scores livrés et leur appariement par sujet et visite avant toute comparaison avec le score M1. Si cet inventaire est incomplet, `cog_scores_mat.json` enregistre le diagnostic et les exports démographiques restent disponibles. Utiliser l’environnement scientifique du cluster et installer les dépendances de `env/requirements-1000brains.txt`, dont `openpyxl` pour le classeur.

```bash
cd /home/aurelio/Server/Aurelio_Stage_M2
source "$HOME/.venvs/lnca-1000brains/bin/activate"
META_RUN=$(date -u +%Y%m%dT%H%M%SZ)
python3 script/metadonnees.py --config config/suivi.json --run "$META_RUN"
python3 script/cognition.py --config config/cognition.json --run "$META_RUN"
```

Le premier import écrit `demographie.csv` (une ligne par sujet et visite), `suivi.csv`, `conflits.csv`, `resume.json` et `cog_scores_mat.json` dans `output/01_population/01_metadonnees/$META_RUN/data/`. Les champs démographiques contradictoires sont laissés vides dans la table et signalés dans `conflits.csv`. Le second écrit `scores.csv` et l’ACP dans `output/01_population/02_scores_cognitifs/$META_RUN/`. Garder les tables par participant sur le cluster, hors Git.

### Lancer et assembler les analyses d’imagerie

1. Préparer un JSON par enregistrement en renseignant les chemins réels. Compléter `config/cohorte.csv` (une ligne par **sujet × session × version**). La colonne `config` est relative au dossier du CSV ou absolue ; `run` est l’identifiant de calcul pour cette ligne. Un doublon est signalé. Commencer par une seule visite en V4. Le lanceur de cohorte ne soumet aucun job automatiquement.
2. Exécuter `individu --row N` pour chaque ligne via les moyens de soumission du laboratoire. Le numéro de ligne commence à 1 après l’en-tête. Cette commande reprend les quatre calculs validés du pilote. Elle réutilise une étape déjà terminée si sa configuration, son code et ses entrées concordent.
3. Assembler les mesures. Pour les figures descriptives selon l’âge, lancer `groupes`. Les options `--reseaux` et `--mc` ajoutent les graphes sur les matrices moyennes (<55 et ≥55 ans), les correspondances entre communautés et les MC moyennes par tertile d’âge. Une FCD individuelle proche de l’âge médian de chaque tertile est illustrée. Les âges et les nombres de participants sont recalculés sur les enregistrements disponibles.
4. Utiliser `--demographie` et `--scores` pour joindre les tables par **sujet et visite**. Le lanceur normalise les identifiants (par exemple `1` et `sub-0001`). L’âge et l’éducation démographiques sont prioritaires quand ils sont renseignés. Les valeurs du classeur cognitif complètent les cellules vides. Les écarts entre sources et les appariements sont comptabilisés dans `jointure_metadonnees.json`. Remplir les colonnes `age,education,cognition` du manifeste seulement pour des valeurs saisies manuellement et vérifiées.
5. Lancer `modeles` sur les sujets complets. Les sorties incluent R² d’ajustement, R² hors pli, prédictions par sujet, ΔR² de l’imagerie au-delà de l’âge et de l’éducation, bootstrap avec réajustement et permutation des résidus du modèle réduit. Examens préalables : distributions et valeurs manquantes, mesures comparables, scores définis et échangeabilité des sujets. Les cinq plis sont constitués par sujet dans une seule session. Lorsque le score cognitif est calculé sur l’ensemble de la table, son échelle de référence est commune ; une prédiction prospective devra définir la référence cognitive hors échantillon ou la réapprendre dans chaque pli.

```bash
python3 script/cohorte.py individu --table config/cohorte.csv --row 1
RUN=$(date -u +%Y%m%dT%H%M%SZ)
DEMO="output/01_population/01_metadonnees/$META_RUN/data/demographie.csv"
SCORES="output/01_population/02_scores_cognitifs/$META_RUN/data/scores.csv"
python3 script/cohorte.py assembler --table config/cohorte.csv --run "$RUN" --demographie "$DEMO" --scores "$SCORES"
python3 script/cohorte.py groupes --table config/cohorte.csv --run "$RUN" --demographie "$DEMO" --scores "$SCORES" --reseaux --mc
python3 script/cohorte.py modeles --table config/cohorte.csv --run "$RUN" --demographie "$DEMO" --scores "$SCORES"
```

| Diapositives M1 | Où retrouver l’analyse |
|---|---|
| 28–32 · population et cognition | Manifeste, distributions d’âge/cognition, ACP de `cognition.py`. |
| 33–45 · structure des réseaux | `statique`, `graphes`, `distance`, `rich-club`, puis `groupes --reseaux`. La correspondance des modules est exportée en table par région et en matrice de comptes. |
| 46–51 · dynamique | `m1`, `modules-mc`, `groupes --mc`. Le sujet pilote et les moyennes de groupe ont des sorties distinctes. |
| 52–60 · projections, tendances et modèles | ACP descriptive, `groupes`, `modeles`. Les options exploratoires M1 (Isomap, interactions spécifiques et figures de résidus) figurent dans l’inventaire historique. Leur reproduction graphique à l’identique reste à adapter à la table définitive. |

## 5. Récupérer les figures et consigner le lancement

Les étapes individuelles écrivent dans `output/<bloc>/<étape>/<RUN>_<sujet>_<session>_<version>/`. Les étapes de population utilisent `output/01_population/<étape>/<RUN>/`. Les manifestes et les journaux permettent de retrouver les entrées et la configuration de chaque calcul.

Pour récupérer l’archive complète, copier l’identifiant UTC affiché par le lanceur, puis exécuter dans PowerShell sur le PC :

```powershell
$run = Read-Host "Identifiant UTC du lancement"
$dest = Join-Path 'C:\Users\aure6\Output_M2Stage' $run
New-Item -ItemType Directory -Force -Path $dest | Out-Null
scp "cluster-lnca:/home/aurelio/Server/Aurelio_Stage_M2/output/04_exports/$run/*.zip" "$dest"
```

Consigner le RUN, le sujet, la session et la décision prise après inspection dans le cahier de laboratoire.

## 6. Suivi longitudinal · réservé pour plus tard

`config/suivi.json` utilise les colonnes de présence FC et SC des visites 1 et 2. Le fichier `suivi.csv` récapitule les visites déclarées. Le champ `colonnes_source` de `resume.json` permet de repérer les éventuelles dates de scan ou un délai. Après vérification du sens de ces champs, renseigner `date_scan` dans les visites et `format_date` si le format diffère de `AAAA-MM-JJ`. Pour un délai fourni directement, renseigner `colonne_intervalle` et `unite_intervalle` (`jours`, `mois` ou `ans`). Les délais sans source restent vides.

Une fois le suivi engagé, vérifier l’appariement avec les fichiers présents et appliquer le même profil aux deux visites. Pour chaque mesure commune, calculer **scan 2 − scan 1** et, lorsque le délai est connu, exprimer ce changement par année. Conserver les sujets sans seconde visite dans le décompte de la cohorte initiale.

## Vérifications techniques

Le test local sur séries synthétiques de 100 régions a vérifié la FC et la FCD contre les corrélations directes, les paires de vitesses, le réalignement d’une SC permutée, un masque avec une interruption et la préservation de la FC par le substitut de phase sur une série continue. Le pilote, la comparaison V1/V4, les figures et un ajustement miniature MFM/tMFM avec évaluation ont été exercés. Les calculs de rich club et de modules de MC ont été testés à taille réduite. La comparaison aux sorties MATLAB et la mesure des ressources du nouveau lancement sur le cluster restent à documenter.

