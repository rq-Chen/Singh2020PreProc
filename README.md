# Preprocessing Code for MINDy

Author: Matthew Singh.

Modified by Ruiqi Chen (2023).

`SuperMessy/` contains all codes. We tried to built the dependency tree and move each package into specific folder. We put all dependencies into the `[PackageName]/utilities` directory. Apart from those, you also need to download `fieldtrip` and specify the fieldtrip directory at the beginning of `MyStartServerDT.m`. You also need `FSL` if DVARS files are not available.

Current availabel HCP data on the WUSTL NIL cluster:

```bash
/net/10.20.145.164/HCPpackages02/unzip/1200subject
/net/10.20.145.162/HCPpackages03/unzip/1200subject
/net/10.20.145.162/HCPpackages04/unzip/1200subject
```

There should be six packages, each with 185-186 subjects. But currently only three are available (as for 2023/06/06).

## For Resting State Only (Singh et. al., 2020)

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

### `MyStartServerDT`

Main function. Documentation is OK.

#### Inputs

- `Subjects`: subject ID, numerical array or cell array of chars
- `nGSR`: version of (Siegel et. al., 2017) preprocessing, set as 3.
- `copyBase`: set as `'y'`
- `Atlas`: set as `'y'` (Yeo network).

#### Pipeline

- Copy the already processed subject folders from `out_base` to `out_dir` (`out_dir` is the final output destination)???
- 

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





## Questions

I am trying to put everything needed (except the HCP data) in a single folder in case some files in `/scratch2` will be deleted one day. Therefore I am trying to understand all paths you include.

For `MyStartServerDT`:

1. What does the `if strcmpi(copyBase(1),'y')` part do? Seems like you are moving some subject folders from `'/scratch2/Singh/HCP/Out/'` to `'/scratch2/Singh/HCP/GSR3/Out/'`. Is there anything special about these subjects? Can I remove this part?
2. You add `/scratch2/Singh/HCP` to path but I only found one `ThalDense.mat` there. Did you use this file? If it's needed, can I just move it to somewhere else and remove the call to `addpath('/scratch2/Singh/HCP')`?
3. You add `/scratch1` to path too but I did not find any matlab files there. Can I remove the addpath?
4. 