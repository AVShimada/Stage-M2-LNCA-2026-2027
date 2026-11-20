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

- Exported the recovered visit 1 cognitive data to CSV and transferred the tables to the cluster.
- Computed the V4 FCD with 30-volume windows (267 windows) and extracted the 10 empirical targets for fitting the mechanistic twin.
- Recomputed the V4 FCD with 10-volume windows and a one-volume step (287 windows) to match the M1 settings.
- Plotted the SIFT2 structural connectivity and the V1 and V4 static FC matrices. Used a shared colour scale to compare FC values.
- Compared V1 and V4 FCD on the same scan with identical window settings. The median similarity between non-overlapping windows was 0.1985 for V1 and 0.1440 for V4.