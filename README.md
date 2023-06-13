# Preprocessing Code for MINDy

Author: Matthew Singh.

Modified and documented by Ruiqi Chen (2023).

This repo contains the fMRI preprocessing codes in the style of (Singh et. al., 2020). There are two versions, one for resting state data only, and the other for resting state and task data simultaneously. Currently we have only cleaned up the one for resting state data and put them into `Singh2020PreProc`. The other scripts are dumped into the folder `SuperMessy`.

HCP data on the WUSTL NIL cluster:

```bash
/net/10.27.136.121/hcpdb/packages/unzip/HCP_1200
```

HCP data on the CHPC cluster:

```bash
/ceph/hcpdb/packages/unzip/HCP_1200
```

They contain all 1113 subjects with neuroimaging data.

## Notes on DVARS

The original (Singh et. al., 2020) paper uses FSL and a bash script (`Singh2020PreProc/utilities/dvars_nichols.sh`) from Dr. Thomas Nichols to compute "standardized DVARS". Seems like this script came from [here](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/fsl/dvars.sh) as given by [Dr. Nichols' 2013 paper](https://arxiv.org/abs/1704.01469).

However, we now use a package developed by Dr. Nichols and Soroosh Afyouni in 2017 called [DVARS](https://github.com/asoroosh/DVARS) associated with their [NeuroImage paper](https://doi.org/10.1016/j.neuroimage.2017.12.098) to calculate the "relative DVARS". For whatever reason, the results from the two methods are different by ~0.1% even after scaling (note that scaling DVARS will not influence our results). We still decided to use this one since it's newer and listed on [Dr. Nichols' website](http://www.nisox.org/Software/DSE/). It also does not require FSL (though it requires MATLAB's statistical toolbox). And most importantly, it has much fewer file IO, making it much faster than the old one on our server (10 min per run -> 100s per run).

However, if you want to use the identical pipeline as Singh2020, you can set `use_old_fslDVARS = true` in `calc_DVAR_mod.m`.

## For Resting State Only (Singh et. al., 2020)

The directory `Singh2020PreProc` contains the main scripts. You'll need to add this folder to MATLAB path.

The subdirectory `utilities` contains the external package `DVARS` and some atlas files.

Apart from those, you need to download `fieldtrip` and specify the fieldtrip directory at the beginning of `MyStartServerDT.m`.

If you want to compute DVARS with the old bash scripts instead of the MATLAB package, you'll also need to set up FSL.

### Preprocessing Pipeline

- HCP protocol: ICA-FIX
- (Siegel et. al., 2017) protocol:
    - Detrending
    - Motion scrubbing:
        - Framewise Displacement (FD)
        - Temporal derivative of variation (DVARS)
    - Respiratory artifact removal
    - Additional motion removal (three versions - the most complete one was used in (Singh et. al., 2020)):
        - 12 HCP motion regressors & their derivatives
        - (on top of above) CompCor: PC of white matter and cerebrospinal fluid signals
        - (on top of above) Global signal regression (GSR): mean signals from white matter, cerebrospinal fluid and grey matter
- Parcellation
- Extra preprocessing by MINDy functions (not included in this repo).

## For Both Resting State & Tasks

### Preprocessing Pipeline

This version does not use ICA-FIX and DVARS-based cut, since those are usually only used for resting state data. Using this version, the resting state and task data will be preprocessed in a same way.

### `MyStartServerDT_HCP_Simple_00`

Main function. Just a wrapper. Calls `MyStartServerDT_HCP_All`.

### `MyStartServerDT_HCP_All`

- Add some paths
- Find or generate and then save `RestOut`.
    - `RestOut` is computed by `MyStartServerDT_HCP_All_restpart`.
- Find or generate and then save `TaskOut` and merge all `TaskOut` into `FullTaskOut`.
    - `TaskOut` is computed by `MyStartServerDT_HCP_All_taskpart`.
- Call `Myc_fcprocess_HCP_All` with `RestOut` and `FullTaskOut`.