function[Out]=MyStartServerDT_HCP_All_Retest_taskpart(Subjects,nGSR,copyBase,Task)

%This Encapsulates the Siegels Resting state Functional Connectivity analysis 
%for DMCC or HCP Subjects

%Basic rundown:
%1.)Build Output file Stucture
%2.)Find or Calculate Both DVars and FDs
%3.)Filter FDs and DVars based on the filter setting below
%4.)Get the WM CSF GM and CompCor regressors
%5.)Setup the Parcellation information
%6.)Either Run the ptseries FC or dtseries FC

%You Need to have the following:
% aparc+aseg.nii.gz for each subject
% _hp2000_clean.nii.gz for every series you plan on running
% _Atlas_hp2000_clean.dtseries.nii for every series
% Movement_Regressors.txt for each series
% Movement_Regressors_dt.txt for each series

%You can run multiple subjects at once through this pipeline
%but do not mix DMCC and HCP subjects because they have different
%parameters
%Make sure that before you run the pipeline all the parameters are setup
%correctly especially the switches 

%the find_FD's and DVAR functions will try to locate the FD's and dvars in
%both the input and the output if the files cannot be found then it will
%generate them and place them in the output folder

%WARNING: 
%Generating the DVARs in this program will take a very ling time
%(about ten minutes per scan) so if you are doing a lot of subjects that
%dont have DVARs you may want to use the cluster to generate DVARs

%OTHER WARNING: 
%if you generate the FD's in this pipeline it assume the movement
%regressors are in degrees and does a conversion this is only true for HCP
%subjects

%The gordon Parcellation used in Siegels original Code was reordered to
%weave together the right and the left hemispheres, this is not usually done
%in parcels, so a function was created that will reorder parcels based on a csv
%look at the function 'reorganize.m' for more info on how to reorder the
%parcel to your liking

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

if ~iscell(Subjects)
    Subjects=num2cell(Subjects);
end
if ~ischar(Subjects{1})
Subjlist=cellfun(@(xx)(num2str(xx)),Subjects,'UniformOutput',0);
else
    Subjlist=Subjects;
end

if strcmpi(Task(1),'W')||strcmpi(Task(1),'N')
    TaskName='WM';
elseif strcmpi(Task(1),'M')
    TaskName='MOTOR';
else
    TaskName=Task;
end
    


RrName={'LR','RL'};
BbName=strcat('tfMRI_',TaskName,'_');
VarNames={'FD_FILT.txt','FD.txt','DVARS.txt','DVARS_FILT.txt','WM_CSF_GM_regs.txt','WM_CSF_CompCor_regs.txt'};
out_base=strcat('/scratch2/Singh/HCP_All_Retest/Out/');
out_dir=strcat('/scratch2/Singh/HCP_All_Retest/GSR',num2str(nGSR),'/Out');

if strcmpi(copyBase(1),'y')
    for ii=1:numel(Subjlist)
        for jj=1:numel(RrName)  %% Scan Sessions
            tExtend=strcat('/',Subjlist{ii},'/Results/',BbName,RrName{jj});

            tDestin=strcat(out_dir,'/',tExtend);
            %tSource=strcat(out_base,tExtend);
            mkdir(tDestin);
            for kk=1:numel(VarNames)
                if kk<5
                    cc=strcat('/',Subjlist{ii},'_',BbName,RrName{jj},'_',VarNames{kk});
                else
                    cc=VarNames{kk};
                end
                copyfile(strcat(out_dir,tExtend,cc),strcat(out_dir,tExtend,cc));
            end
        end
    end
end


%clear
%in_dir='/scratch2/Singh/HCP/';
%in_dir='/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/MINIMALLY_PREPROCESSED';

%in_dir  = '/data/hcp-zfs/OpenAccess/900subject/'; % the Directoy where your subject data lives
in_dir='/net/zfs-HCPpackages04/HCPpackages/unzip/HCP_Retest/';
%out_dir = '/scratch1/MitchJeffers/restingState/out/HCP'% Where the results will be output
%Subjlist = {'132017', '155938','568963','594156' }%list of subjects to be ran through the pipeline
%Subjlist={'100307'};
%in_dir='C:\Users\Matthew\Desktop\HCP';
%out_dir='C:\Users\Matthew\Desktop\HCP\Siegel';


%% add external tools

%ToolPath='C:\Users\Matthew\Desktop\Control\Resting Processing\';
%toolNames={'FSLNets','fieldtrip','fieldtrip\utilities','fieldtrip\external\freesurfer','fieldtrip\fileio','fieldtrip\fileio\special'};

%for i=1:numel(toolNames)
%    addpath(strcat(ToolPath,toolNames{i}));
%end
addpath('/scratch2/Singh/HCP')
%addpath(genpath('/scratch1/MitchJeffers/restingState/tests/TEST'))
addpath('/scratch1/MitchJeffers/restingState/FSLNets')          % for function nets_netmats
addpath('/scratch1/MitchJeffers/restingState/fieldtrip')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/utilities')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/external/freesurfer')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio/special') % (some fieldtrip stuff used by Siegel..not mentioned, though)
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio')
addpath('/scratch1')


%Using HCP data or DMCC data?
Data.Type ='HCP';%HCP or DMCC
tseries = RestNameSetup_Task(Task);%sets up the tseries name used in either DMCC or HCP
switches = SwitchSetup(Data.Type);%sets up the parameters CHECK BEFORE RUNNING
switches.GSR=nGSR;
switches.FDcut=.9;
switches.DVARcut = 1.5;
switches.minframes = 5;

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
    BuildFileStructure(out_dir, Subjlist{i}, tseries)
    
%        for jj=1:numel(RrName)  %% Scan Sessions
%            tExtend=strcat('/',Subjlist{i},'/Results/',BbName,RrName{jj});
%            BadFilesAll(tExtend,100);
%        end

    
    
    %have the DVARs or FD's been Calculated?
    %WARNING If you try to make DVARs on CCPLINUX1 it will take a really 
    %long time.(like 10 mins per series) you may be better off loading them on to the cluster if you
    %have a lot of subjects to do
    find_FD(in_dir, out_dir, Subjlist{i},tseries);
    find_DVAR(in_dir, out_dir, Subjlist{i}, tseries);

    %Filter the FD's and DVARS using the predefined 'Filter' to remove
    %respritory
    %to do fix where to look for HCP subjects
    disp('Filtering FDs and DVARS')
    filter_FD_DVARS(Filter, Subjlist{i}, in_dir, out_dir,tseries)

    
    %Get the regressors for all the functioal connectivity
    disp('Running HCPregressorCompCor')
    getHCPregressorsCompCor_final_NoFIX(Subjlist{i}, in_dir, out_dir,tseries)
end
 

%If you are running the Ptseries choose which parcellation you want to use
% Parc.dir = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/HCP-MMP/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'
% Parc.NP = 360;%Number of Parcels
% Parc.name = 'MMP360';%Name of the parcel will be used to store the ptseries
% Parc.ordered = 0;%is the parcel already ordered Correctly ie. Community then Hemisphere
% Parc.key = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/HCP-MMP/Glasser_et_al_2016_HCP_MMP1.0_RVVG/MMP360ParcelsKey.csv'

%ParDir='C:\Users\Matthew\Desktop\HCP\ATLASES\';
%Parc.dir=strcat(ParDir,'gordon\gordon_parcels\Parcels\Parcels_LR.dlabel.nii');


Parc.dir = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/gordon/gordon_parcels/Parcels/Parcels_LR.dlabel.nii';%where is the parcellation located
Parc.NP = 333; % number of parcels in Gordon parcellation, change accordingly for other parcellations
Parc.name = 'Gordon333'; %what should we call the output of the ptseries .mats
Parc.ordered = 0; % is the parcel ordered on a community level
%Parc.key=strcat(ParcDir,'gordon\gordon_parcels\Parcels\Gordon333Parcelskey.csv');
%Parc.key = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/gordon/gordon_parcels/Parcels/Gordon333ParcelsKey.csv'%path to the key of how the parcels should be reordered
Parc.key = '/scratch1/MitchJeffers/restingState/Gordon333SRI_ParcelsKey.csv'; %#ok<STRNU>

%Parc.dir = '/data/nil-bluearc/corbetta/PP_SCRIPTS/Parcellation/GLParcels/reordered/GLParcels_324_reordered.32k.dlabel.nii'
%Parc.NP = 324
%Parc.name = 'Gordon324'
%Parc.ordered = 1

 disp('running dtSeries')%run a dense timeseries
% Myc_fcprocess_HCP_dtseriesTask(switches, Subjlist, in_dir, out_dir, tseries,Task)
%disp('running ptSeries')%run parcellated time series
%c_fcprocess_HCP_ptseries(switches, Subjlist, in_dir, out_dir, tseries, Parc)
%d112817ParcellateServer(str2double(Subjlist),nGSR,[],Atlas);
switches.FIX=false;
switches.DVARcut=false;
Out={switches, in_dir, out_dir, tseries,Task};

end
