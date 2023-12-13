# Preprocessing Code for MINDy Resting State Model

Author: Matthew Singh.

Modified and documented by Ruiqi Chen (2023).

This repo contains the fMRI preprocessing codes in the style of (Singh et. al., 2020). Here is the pipeline (for details see (Singh et al., 2020)):

- HCP protocol: ICA-FIX
- (Siegel et. al., 2017) protocol:
    - Detrending
    - Motion scrubbing:
        - Framewise Displacement (FD): >0.2mm
        - Temporal derivative of variation (DVARS): >1.05*median.
            - Note: DVARS was calculated using the "minimally preprocessed" data from HCP WITHOUT ICA-FIX!
        - FD and DVARS were filtered for respiratory artifact before use
        - Note: the last frame was discarded.
    - Additional nuisance regression (three versions - the most complete one was used in (Singh et. al., 2020)):
        - no extra removal
        - (on top of above) CompCor: the top five PCs of white matter and cerebrospinal fluid signals
        - (on top of above) Global signal regression (GSR): mean signal from gray matter
- Parcellation according to all levels of Schaefer atlas (Schaefer et al., 2018):
    - Note: the first frame was discarded.
- Extra preprocessing by MINDy functions (not included in this repo).

## Notes on DVARS

The original (Singh et al., 2020) paper uses FSL and a bash script (`Singh2020PreProc/utilities/dvars_nichols.sh`) from Thomas Nichols to compute "standardized DVARS". Seems like this script came from [here](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/fsl/dvars.sh) as given by [Thomas Nichols' 2013 paper](https://arxiv.org/abs/1704.01469).

However, we now use a package developed by Nichols and Afyouni in 2017 called [DVARS](https://github.com/asoroosh/DVARS) associated with their [NeuroImage paper](https://doi.org/10.1016/j.neuroimage.2017.12.098) to calculate the "relative DVARS". For whatever reason, the results from the two methods are different by ~0.1% even after scaling (note that scaling DVARS will not influence our results). We still decided to use this one since it's newer and listed on [Thomas Nichols' website](http://www.nisox.org/Software/DSE/). It also does not require FSL (though it requires MATLAB's statistical toolbox). And most importantly, it has much fewer file IO, making it much faster than the old one on our server (10 min per run -> 100s per run).

However, if you want to use the identical pipeline as Singh2020, you can set `use_old_fslDVARS = true` in `calc_DVAR_mod.m`.

## Directory Structure

The directory `Singh2020PreProc` contains the main scripts. You'll need to add this folder to MATLAB path. Alternatively, you can call the wrapper `MyWrapper.m` from the root which will identify the unprocessed data folders, manage the paths and run the codes.

The subdirectory `Singh2020PreProc/utilities` contains the external package `DVARS` and some atlas files.

Apart from those, you need to download `fieldtrip` and specify the fieldtrip directory at the beginning of `Singh2020PreProc/MyStartServerDT.m`.

If you want to compute DVARS with the old bash scripts instead of the MATLAB package, you'll also need to set up FSL.

## Data

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

Singh, M. F., Braver, T. S., Cole, M. W., & Ching, S. (2020). Estimation and validation of individualized dynamic brain models with resting state fMRI. NeuroImage, 221, 117046. https://doi.org/10.1016/j.neuroimage.2020.117046

Siegel, J. S., Mitra, A., Laumann, T. O., Seitzman, B. A., Raichle, M., Corbetta, M., & Snyder, A. Z. (2017). Data Quality Influences Observed Links Between Functional Connectivity and Behavior. Cerebral Cortex, 27(9), 4492–4502. https://doi.org/10.1093/cercor/bhw253

Schaefer, A., Kong, R., Gordon, E. M., Laumann, T. O., Zuo, X.-N., Holmes, A. J., Eickhoff, S. B., & Yeo, B. T. T. (2018). Local-Global Parcellation of the Human Cerebral Cortex from Intrinsic Functional Connectivity MRI. Cerebral Cortex, 28(9), 3095–3114. https://doi.org/10.1093/CERCOR/BHX179