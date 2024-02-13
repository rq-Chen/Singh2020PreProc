# Preprocessing Code for MINDy HCP rfMRI Model

Preprocessing pipeline as used in (Singh et al., 2020) and (Chen et al., 2024).

Author: Joshua Siegel and Matthew Singh.

Modified and documented by Ruiqi Chen.

This repo contains the preprocessing scripts for HCP rfMRI data in the style of (Singh et al., 2020). Here is the pipeline (for details see (Singh et al., 2020)):

- HCP protocol: FIX-cleaned data (with FIX-ICA and motion correction)
- (Siegel et al., 2017) protocol:
    - Detrending
    - Motion scrubbing (linearly interpolating (but not extrapolating) high-motion frames):
        - Framewise Displacement (FD): >0.2mm
        - Temporal derivative of variation (DVARS): >1.05*median.
        - FD and DVARS were filtered for respiratory artifact before use
    - Additional nuisance regression (three versions - the most complete one was used in (Singh et al., 2020)):
        - (`GSR0`) just remove parcel mean;
        - (`GSR2`) (on top of above) CompCor: the top five PCs of white matter and cerebrospinal fluid signals;
        - (`GSR3`) (on top of above) Global signal regression (GSR): mean signal from gray matter.
- Averaging within each parcel accroding to all levels of Schaefer atlas (Schaefer et al., 2018):
- The first frame was discarded.
- Extra preprocessing by MINDy functions (not included in this repo).

Also note that although the scripts can be used to preprocessed both "minimally preprocessed" and FIX-ICA versions of data, the DVARS is always calculated using the "minimally preprocessed" (no FIX-ICA) data. This is more conservative because more high-motion frames will be identified and interpolated when using this more noisy version of data.

The entry point is `Singh2020PreProc/MyStartServerDT`. See the comments in that file for details. You need to add that folder to MATLAB path.

## Notes on DVARS

The original (Singh et al., 2020) paper uses FSL and a bash script (`Singh2020PreProc/utilities/dvars_nichols.sh`) from Thomas Nichols to compute "standardized DVARS". Seems like this script came from [here](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/fsl/dvars.sh) as given by [Thomas Nichols' 2013 paper](https://arxiv.org/abs/1704.01469).

However, we now use a package developed by Nichols and Afyouni in 2017 called [DVARS](https://github.com/asoroosh/DVARS) associated with their [NeuroImage paper](https://doi.org/10.1016/j.neuroimage.2017.12.098) to calculate the "relative DVARS". For whatever reason, the results from the two methods are different by ~0.1% even after scaling (note that scaling DVARS will not influence our results). We still decided to use this one since it's newer and listed on [Thomas Nichols' website](http://www.nisox.org/Software/DSE/). It also does not require FSL (though it requires MATLAB's statistical toolbox). And most importantly, it has much fewer file IO, making it much faster than the old one on our server (10 min per run -> 100s per run).

However, if you want to use the identical pipeline as Singh2020, you can set `use_old_fslDVARS = true` in `calc_DVAR_mod.m`.

## Directory Structure

The directory `Singh2020PreProc` contains the main scripts. You'll need to add this folder to MATLAB path. Alternatively, you can call the wrapper `MyWrapper.m` from the root which will identify the unprocessed data folders, manage the paths and run the codes.

The subdirectory `Singh2020PreProc/utilities` contains the external package `DVARS` and some atlas files.

Apart from those, you need to download `fieldtrip` and specify the fieldtrip directory at the beginning of `Singh2020PreProc/MyStartServerDT.m`.

If you want to compute DVARS with the old bash scripts instead of the MATLAB package, you'll also need to set up FSL.

### I/O

The default input and output folders are specified in `MyWrapper.m`. The input folder should be the HCP folder (containing one subfolder for each participant) and the output folder will contain one MAT file for each participant and each level of parcellation, entitled like `sub[HCP_ID]Y[xx].mat` where `[xx]` is the number of parcels divided by 100, e.g., `sub100206Y01.mat`. Each file contains a single structure `X` with the following fields that worth mentioning:

- `Dat`: The data. `(1, nRuns)` cell array with each cell being a `(nParcels + 19, 1199)` matrix. `nRuns` is usually four (`{'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}`) but can be smaller for participants with missing data. The parcels are in the order of Schaefer atlas (Schaefer et al., 2018), followed by the 19 FreeSurfer subcortical regions. The timepoint corresponds to TR `2:1200` in the original data.
- `GSR`: Level of nuisance regression. Default is 3.
    - `GSR = 0`: `regs = ones(1200, 1)`;
    - `GSR = 1` (not mentioned above): `regs = [mean_white_matter mean_CSF mean_gray_matter ones(1200, 1)]`;
    - `GSR = 2`: `regs = [Top5PC_white_matter Top5PC_CSF ones(1200, 1)]`;
    - `GSR = 3`: `regs = [Top5PC_white_matter Top5PC_CSF mean_gray_matter ones(1200, 1)]`;
- `QC`: Quality control statistics. A structure with the following fields:
    - `frames`: `(1, nRuns)` matrix, number of motion-free frames in each run.
    - `run`: `(1, nRuns)` structure array, with the following fields:
        - `tmask`: `(1, 1200)` binary vector, whether each frame should be excluded based on movement.
        - `FD`: `(1, 1200)` vector, the framewise displacement.
        - `DV`: `(1, 1200)` vector, the DVARS.
    - `regs`: `(1200, nReg)` matrix, the regressors used for nuisance regression. `nReg` depends on `GSR`.

## Data (WUSTL only)

HCP data on the WUSTL NIL cluster:

```bash
/net/10.27.136.121/hcpdb/packages/unzip/HCP_1200
```

HCP data on the CHPC cluster:

```bash
/ceph/hcpdb/packages/unzip/HCP_1200
```

They contain all 1113 subjects with neuroimaging data.

## References

Chen, R., Singh, M., Braver, T. S., & Ching, S. (2024). Dynamical models reveal anatomically reliable attractor landscapes embedded in resting state brain networks. bioRxiv. [https://doi.org/10.1101/2024.01.15.575745](https://doi.org/10.1101/2024.01.15.575745)

Singh, M. F., Braver, T. S., Cole, M. W., & Ching, S. (2020). Estimation and validation of individualized dynamic brain models with resting state fMRI. NeuroImage, 221, 117046. [https://doi.org/10.1016/j.neuroimage.2020.117046](https://doi.org/10.1016/j.neuroimage.2020.117046)

Siegel, J. S., Mitra, A., Laumann, T. O., Seitzman, B. A., Raichle, M., Corbetta, M., & Snyder, A. Z. (2017). Data Quality Influences Observed Links Between Functional Connectivity and Behavior. Cerebral Cortex, 27(9), 4492–4502. [https://doi.org/10.1093/cercor/bhw253](https://doi.org/10.1093/cercor/bhw253)

Schaefer, A., Kong, R., Gordon, E. M., Laumann, T. O., Zuo, X.-N., Holmes, A. J., Eickhoff, S. B., & Yeo, B. T. T. (2018). Local-Global Parcellation of the Human Cerebral Cortex from Intrinsic Functional Connectivity MRI. Cerebral Cortex, 28(9), 3095–3114. [https://doi.org/10.1093/CERCOR/BHX179](https://doi.org/10.1093/CERCOR/BHX179)
