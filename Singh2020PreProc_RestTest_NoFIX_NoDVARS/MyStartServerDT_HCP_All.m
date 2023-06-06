function[]=MyStartServerDT_HCP_All(Subjects,nGSR,copyBase,TaskAll,Atlas)
%% TaskAll should be a cell
%% First Rep is Resting-State

%% NOTE: This version doesn't use ICA-FIX or DVAR-based cuts to ensure Task and Rest match
%% However, DVAR time series and cut thresholds are saved so cuts can be done later.

if ~iscell(TaskAll)
    TaskAll={TaskAll};
end
addpath('/scratch2/Singh/HCP')
%addpath(genpath('/scratch1/MitchJeffers/restingState/tests/TEST'))
addpath('/scratch1/MitchJeffers/restingState/FSLNets')          % for function nets_netmats
addpath('/scratch1/MitchJeffers/restingState/fieldtrip')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/utilities')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/external/freesurfer')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio/special') % (some fieldtrip stuff used by Siegel..not mentioned, though)
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio')
addpath('/scratch1')

OutDir=strcat('/scratch2/Singh/HCP_All/GSR',num2str(nGSR),'/Out/');
disp('Doing Rest')
for iSub=1:numel(Subjects)
    if exist(strcat(OutDir,Subjects{iSub},'/RestOut.mat'), 'file') == 2
        RestExist(iSub)=1;
    else
        RestExist(iSub)=0;
    end
    if mean(RestExist)==1
        load(strcat(OutDir,Subjects{iSub},'/RestOut.mat'))
    else
        RestOut=MyStartServerDT_HCP_All_restpart(Subjects(RestExist==0),nGSR,copyBase);
        for iS=find(RestExist==0)
            save(strcat(OutDir,Subjects{iS},'/RestOut.mat'),'RestOut');
        end
    end
    disp('Doing Task')
    FullTaskOut=cell(1,numel(TaskAll));
    for iTask=1:numel(TaskAll)
        disp([iTask numel(TaskAll)])
        for iSub=1:numel(Subjects)
            if exist(strcat(OutDir,Subjects{iSub},'/TaskOut_',TaskAll{iTask},'.mat'), 'file') == 2
                TaskExist(iSub)=1;
            else
                TaskExist(iSub)=0;
            end
        end
        if mean(TaskExist)==1
            load(strcat(OutDir,Subjects{iSub},'/TaskOut_',TaskAll{iTask},'.mat'))
        else
            TaskOut=MyStartServerDT_HCP_All_taskpart(Subjects(TaskExist==0),nGSR,copyBase,TaskAll{iTask});
            for iS=find(TaskExist==0)
                save(strcat(OutDir,Subjects{iS},'/TaskOut_',TaskAll{iTask},'.mat'),'TaskOut');
            end
        end
        FullTaskOut{iTask}=TaskOut;
    end
end

disp('running dtSeries')%run a dense timeseries

Myc_fcprocess_HCP_All(nGSR,Subjects,RestOut,FullTaskOut,TaskAll,Atlas);
%Myc_fcprocess_HCP_dtseriesTask(switches, Subjlist, in_dir, out_dir, tseries,Task)
%disp('running ptSeries')%run parcellated time series
%c_fcprocess_HCP_ptseries(switches, Subjlist, in_dir, out_dir, tseries, Parc)
%d112817Parce

end
