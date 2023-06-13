function[]=MyStartServerDT(Subjects, nGSR, Atlas, out_dir, in_dir)
%
% Update by Ruiqi Chen 2023/06/06:
%
% Inputs:
%   - Subjects: subject IDs (e.g., 100307), numerical array or cell of char
%   - nGSR: level of (Siegel et. al., 2017) preprocessing. Set as 3 for (Singh
%   et. al., 2020) pipeline.
%   - Atlas: 'y' for Yeo17 (Schaefer) atlas, otherwise Gordon333 (not implemented yet)
%   - out_dir: Output directory
%   - in_dir: HCP directory, with subdirectories named after subject IDs.
%
% Outputs: Will be written to the output directory. Contains:
%   - MATLAB files `sub[SubjID]Y[nParcels/100].mat`: e.g., `sub100206Y01.mat`,
%       `sub149741Y10.mat`. These are the inputs to MINDy functions (parcel-level
%       timeseries `X.Dat` along with quality control (QC) measures `X.QC`).
%   - Subdirectory `Out/`:
%     - Subdirecoty `outputs/`: dense timeseries for each subject
%       - one long-named .mat file (size ~2G) for each subject
%     - Subdirectories `[SubjID]/`: e.g., `100307/`
%       - Subdirectory `Results/`:
%         - Subdirectories for each run: e.g., `rfMRI_REST1_LR/`:
%           - DVARS file
%           - filtered DVARS file
%           - FD file
%           - filtered FD file
%
%
%This Encapsulates the Siegels Resting state Functional Connectivity analysis
%for DMCC or HCP Subjects
%
%Basic rundown:
%1.)Build Output file Stucture
%2.)Find or Calculate Both DVars and FDs
%3.)Filter FDs and DVars based on the filter setting below
%4.)Get the WM CSF GM and CompCor regressors
%5.)Setup the Parcellation information
%6.)Either Run the ptseries FC or dtseries FC
%
%You Need to have the following:
% aparc+aseg.nii.gz for each subject
% _hp2000_clean.nii.gz for every series you plan on running
% _Atlas_hp2000_clean.dtseries.nii for every series
% Movement_Regressors.txt for each series
% Movement_Regressors_dt.txt for each series
%
%You can run multiple subjects at once through this pipeline
%but do not mix DMCC and HCP subjects because they have different
%parameters
%Make sure that before you run the pipeline all the parameters are setup
%correctly especially the switches
%
%the find_FD's and DVAR functions will try to locate the FD's and dvars in
%both the input and the output if the files cannot be found then it will
%generate them and place them in the output folder
%
%WARNING:
%Generating the DVARs in this program will take a very long time
%(about ten minutes per scan) so if you are doing a lot of subjects that
%dont have DVARs you may want to use the cluster to generate DVARs
%
%OTHER WARNING:
%if you generate the FD's in this pipeline it assume the movement
%regressors are in degrees and does a conversion this is only true for HCP
%subjects
%
%The gordon Parcellation used in Siegels original Code was reordered to
%weave together the right and the left hemispheres, this is not usually done
%in parcels, so a function was created that will reorder parcels based on a csv
%look at the function 'reorganize.m' for more info on how to reorder the
%parcel to your liking
%
%Outputs:
%FD's and DVARs
%Filtered FDs and DVARs
%WM_CSF_CompCor regressors
%WM_CSF_GM regressors
%.ptseries (if your using the ptseries)
%and a long named .mat file placed in the toplevel outputs folder, it contains:
%HCP_GL_FC which is the Cross run average FC at each parcel for each
%subject ran
%and QC which contains masking, frame and run exclusion, FD and DVAR info.
%as well as the individual run FC


%% Constants %%

fieldtrip_dir = fullfile(fileparts(mfilename('fullpath')), 'utilities', 'fieldtrip')


%% Inputs %%

if ~iscell(Subjects)
    Subjects=num2cell(Subjects);
end
if ~ischar(Subjects{1})
    Subjlist=cellfun(@(xx)(num2str(xx)),Subjects,'UniformOutput',0);
else
    Subjlist=Subjects;
end

script_dir = fileparts(mfilename('fullpath'));
if nargin < 5
    in_dir = '/net/10.20.145.162/HCPpackages03/unzip/1200subject';
    fprintf('\n\nNo input directory provided. Default to %s\n\n', in_dir);
end
if nargin < 4
    out_dir = fullfile(script_dir, '..', ['GSR' num2str(nGSR)]);
    fprintf('\n\nNo output directory provided. ');
    fprintf('Default to a new folder in the father directory of this script: %s\n\n', out_dir);
end
dtseries_dir = fullfile(out_dir, 'Out');

RrName={'1_LR','1_RL','2_LR','2_RL'};
BbName='rfMRI_REST';
% VarNames={'FD_FILT.txt','FD.txt','DVARS.txt','DVARS_FILT.txt','WM_CSF_GM_regs_hp200_clean.txt','WM_CSF_CompCor_regs_hp200_clean.txt'};
% out_base=strcat('/scratch2/Singh/HCP/Out/');
% dtseries_dir=strcat('/scratch2/Singh/HCP/GSR',num2str(nGSR),'/Out');

% if strcmpi(copyBase(1),'y')
%     for ii=1:numel(Subjlist)
%         for jj=1:numel(RrName)  %% Scan Sessions
%             tExtend=strcat('/',Subjlist{ii},'/Results/',BbName,RrName{jj},'/');

%             tDestin=strcat(dtseries_dir,'/',tExtend);
%             %tSource=strcat(out_base,tExtend);
%             mkdir(tDestin);
%             for kk=1:numel(VarNames)
%                 if kk<5
%                     cc=strcat('/',Subjlist{ii},'_',BbName,RrName{jj},'_',VarNames{kk});
%                 else
%                     cc=VarNames{kk};
%                 end
%                 copyfile(strcat(out_base,tExtend,cc),strcat(dtseries_dir,tExtend,cc));
%             end
%         end
%     end
% end

% in_dir  = '/data/hcp-zfs/OpenAccess/1200subject/'; % the Directory where your subject data lives


%% add external tools

% addpath('/scratch2/Singh/HCP')
% addpath(fullfile(script_dir, 'utilities', 'FSLNets'))          % for function nets_netmats
addpath(fullfile(fieldtrip_dir))
addpath(fullfile(fieldtrip_dir, 'utilities'))
addpath(fullfile(fieldtrip_dir, 'external', 'freesurfer'))
addpath(fullfile(fieldtrip_dir, 'fileio', 'special')) % (some fieldtrip stuff used by Siegel..not mentioned, though)
addpath(fullfile(fieldtrip_dir, 'fileio'))
% addpath('/scratch1')


%Using HCP data or DMCC data?
Data.Type ='HCP';%HCP or DMCC
tseries = RestNameSetup(Data.Type);%sets up the tseries name used in either DMCC or HCP
switches = SwitchSetup(Data.Type);%sets up the parameters CHECK BEFORE RUNNING
switches.GSR=nGSR;


% for [0.06 - 0.14] Hz, removing respiratory frequencies (Siegel) Used for
% FDs and DVARs
Filter = designfilt('bandstopfir', ...
    'FilterOrder', 40, ...
    'CutoffFrequency1', 0.06, ...
    'CutoffFrequency2', 0.14, ...
    'SampleRate', 1/switches.TR);

%build the File Structure to store the result data. and build the Filtered
%FD and DVars
for i = 1:length(Subjlist)

    BuildFileStructure(dtseries_dir, Subjlist{i}, tseries)

    for jj=1:numel(RrName)  %% Scan Sessions
        tExtend=strcat(dtseries_dir,'/',Subjlist{i},'/Results/',BbName,RrName{jj},'/');
        BadFilesAll(tExtend,100);
    end

    t0 = tic;

    %have the DVARs or FD's been Calculated?
    %WARNING If you try to make DVARs on CCPLINUX1 it will take a really
    %long time.(like 10 mins per series) you may be better off loading them on to the cluster if you
    %have a lot of subjects to do
    find_FD(in_dir, dtseries_dir, Subjlist{i},tseries);
    find_DVAR(in_dir, dtseries_dir, Subjlist{i}, tseries);

    %Filter the FD's and DVARS using the predefined 'Filter' to remove
    %respritory
    %to do fix where to look for HCP subjects
    disp('Filtering FDs and DVARS')
    filter_FD_DVARS(Filter, Subjlist{i}, in_dir, dtseries_dir, tseries)

    disp(['Elapsed time for generating filtered FD & DVARS for subject ' Subjlist{i} ' :'])
    toc(t0)
    t0 = tic;

    %Get the regressors for all the functioal connectivity
    disp('Running HCPregressorCompCor')
    getHCPregressorsCompCor_final(Subjlist{i}, in_dir, dtseries_dir, tseries)

    disp(['Elapsed time for generating CompCor regressors for subject ' Subjlist{i} ' :'])
    toc(t0)
end


%If you are running the Ptseries choose which parcellation you want to use
% Parc.dir = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/HCP-MMP/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'
% Parc.NP = 360;%Number of Parcels
% Parc.name = 'MMP360';%Name of the parcel will be used to store the ptseries
% Parc.ordered = 0;%is the parcel already ordered Correctly ie. Community then Hemisphere
% Parc.key = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/HCP-MMP/Glasser_et_al_2016_HCP_MMP1.0_RVVG/MMP360ParcelsKey.csv'

%ParDir='C:\Users\Matthew\Desktop\HCP\ATLASES\';
%Parc.dir=strcat(ParDir,'gordon\gordon_parcels\Parcels\Parcels_LR.dlabel.nii');


% Parc.dir = fullfile(script_dir, 'utilities', 'Parcels_LR.dlabel.nii');%where is the parcellation located
% Parc.NP = 333; % number of parcels in Gordon parcellation, change accordingly for other parcellations
% Parc.name = 'Gordon333'; %what should we call the output of the ptseries .mats
% Parc.ordered = 0; % is the parcel ordered on a community level
% %Parc.key=strcat(ParcDir,'gordon\gordon_parcels\Parcels\Gordon333Parcelskey.csv');
% %Parc.key = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/gordon/gordon_parcels/Parcels/Gordon333ParcelsKey.csv'%path to the key of how the parcels should be reordered
% Parc.key = fullfile(script_dir, 'utilities', 'Gordon333SRI_ParcelsKey.csv'); %#ok<STRNU>

%Parc.dir = '/data/nil-bluearc/corbetta/PP_SCRIPTS/Parcellation/GLParcels/reordered/GLParcels_324_reordered.32k.dlabel.nii'
%Parc.NP = 324
%Parc.name = 'Gordon324'
%Parc.ordered = 1

t0 = tic;
disp('running dtSeries')%run a dense timeseries
Myc_fcprocess_HCP_dtseries(switches, Subjlist, in_dir, dtseries_dir, tseries)
%disp('running ptSeries')%run parcellated time series
%c_fcprocess_HCP_ptseries(switches, Subjlist, in_dir, dtseries_dir, tseries, Parc)if strcmpi(Atlas(1),'y')
disp('Elapsed time for dense timeseries preprocessing:')
toc(t0)

t0 = tic;
if strcmpi(Atlas(1),'y')
    % if strcmpi(Atlas(2:end),'7')
    %     AllYeoParcellate7net(str2double(Subjlist),nGSR,'y',[],'n');
    % else
        MyAllYeoParcellate(str2double(Subjlist), nGSR, out_dir, [], 'n');
    % end
else
    % d112817ParcellateServer(str2double(Subjlist),nGSR,[],Atlas);
    error('Not implemented!')
end
disp('Elapsed time for parcellation:')
toc(t0)

end
