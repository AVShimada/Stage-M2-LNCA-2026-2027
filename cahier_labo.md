## Monday and Tuesday, 16–17 November 2026

- Read the internship papers on resting-state dynamics, FCD and mechanistic brain twins (Deco et al., 2011; Hansen et al., 2015; Battaglia et al., 2020; Pathak and Battaglia, 2026).
- Attended two meetings and one class.

## Wednesday, 18 November

- Read the cluster’s Slurm documentation.
- Set up the project on the cluster, generated `architecture_projet.txt`, and connected the Git repository to GitHub.

## Thursday, 19 November

- Copied the 1000BRAINS dataset to the cluster and set up the workflow to reproduce the M1 analyses on one subject.
- Selected `sub-0001 / ses-1 / V4` with the Schaefer 100 atlas. Verified 296 BOLD volumes and checked the agreement between the supplied and recalculated static FC.

## Friday, 20 November

- Prepared the visit 1 cognitive CSVs and recovered the original 1000BRAINS metadata files from the cluster (`Demographic_data.csv`, `1000Brains_NP_VarOI_V1_Pseudo.xlsx`, and `cog_scores.mat`). Stored them in `1kBrains/metadonnees/`.
- Computed the V4 FCD with 30-volume windows (267 windows) and extracted the 10 empirical targets for fitting the mechanistic twin.
- Recomputed the V4 FCD with 10-volume windows and a one-volume step (287 windows) to match the M1 settings.
- Plotted the SIFT2 structural connectivity and the V1 and V4 static FC matrices. Used a shared colour scale to compare FC values.
- Compared V1 and V4 FCD on the same scan with identical window settings. The median similarity between non-overlapping windows was 0.1985 for V1 and 0.1440 for V4.

## Tuesday 24 November 2026

- Reorganized the workflow into population analyses, empirical pilot analyses, twin simulations and exports.
- Selected 409 of the 432 participants with a declared second scan. Excluded 23 participants with at least one missing age-appropriate cognitive test. The AKT was required only when applicable.
- Recomputed demographic and cognitive summaries. The cohort included 179 women and 230 men. Mean age at scan 1 was 61.3 years. The median interval estimated from age at each scan was 3.7 years (408 participants with both ages available).
- Ran the imaging workflow on sub-0001, session 1, using all 296 volumes and the Schaefer 100 atlas. Compared V1 and V4, then extracted empirical FC, FCD and the ten dynamic targets.

## Wednesday 25 November 2026

- Completed the pilot graph analyses, including distance, rich club and meta-connectivity modules.
- Calibrated the fixed MFM. Increasing G from 0.05 to 0.75 at sigma = 0.03 reduced mean validation error on the ten dynamic targets from 1.120 to 0.206 across five simulation seeds.
- Examined the static FC discrepancy. Median empirical FC was 0.012 in V1 and 0.388 in V4. After regressing the global signal from both series as a diagnostic, the V1–V4 FC RMSE fell from 0.392 to 0.086. The original denoising steps still need confirmation.

## Thursday 26 November 2026

- Tested time-varying global coupling on the same pilot. The initial modulation produced near-perfect regional correlations. Repeated the sweep at tau = 10 s across ten amplitudes from 0 to 0.05, including the fixed MFM.
- Selected A = 0.00003 on the fitting seeds. Across five new validation seeds, mean dynamic-target error fell from 0.224 to 0.191 and static FC RMSE from 0.382 to 0.283. Both improved on every seed.
- Retained A = 0.00003 as the pilot reference. Median simulated FC remained below the empirical value (0.188 vs 0.388), and the spatial correlation between simulated and empirical FC profiles did not improve.