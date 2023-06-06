function Myc_fcprocess_HCP_All(GSR,Subjlist,RestC,TaskCell,TaskAll,Atlas)

nSubs=numel(Subjlist);
RestMean=cell(1,nSubs);
RestSTD=cell(1,nSubs);
TaskMean=cell(1,nSubs);
TaskSTD=cell(1,nSubs);

excluded = zeros(1,1);

%addpath('/scratch2/Singh/HCP')
%addpath(genpath('/scratch1/MitchJeffers/restingState/tests/TEST'))
%addpath('/scratch1/MitchJeffers/restingState/FSLNets')          % for function nets_netmats
%addpath('/scratch1/MitchJeffers/restingState/fieldtrip')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/utilities')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/external/freesurfer')
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio/special') % (some fieldtrip stuff used by Siegel..not mentioned, though)
addpath('/scratch1/MitchJeffers/restingState/fieldtrip/fileio')
%addpath('/scratch1')

ParcMask=GetAtlasMasks(Atlas,'y');

if strcmpi(Atlas(1),'Y')
    nParcLevels=numel(ParcMask);
else
    nParcLevels=1;
end

%% Do the stuff
for ll = 1:(nSubs) % subjects
    TaskName=cell(1,numel(TaskCell)+1);
    FullQC=cell(1,numel(TaskCell)+1);
    FullDat=cell(nParcLevels,numel(TaskCell)+1);

    for iTask=1:(numel(TaskAll)+1)

        QC = [];
        regs = [];
        runs = zeros(nSubs,1);
        if iTask==1
            switches=RestC{1};source_dir=RestC{2};work_dir=RestC{3};tseries=RestC{4};
        else
            TaskC=TaskCell{iTask-1};
            switches=TaskC{1};source_dir=TaskC{2};work_dir=TaskC{3};tseries=TaskC{4};
        end
        NR = length(tseries);   % number of runs (LR and RL in HCP)

        switches
        %% Load a bunch of settings (You can preset them in the switches struct)

        if ~isfield(switches, 'AveRun');
            switches.AveRun = input('Do you want to concatenate timecourses (0) or calc FC by run and then average (1)?');
        end

        if ~isfield(switches, 'GSR');
            switches.GSR = input('Timecourse Regression? (0 = No, 1 = GSR, 2 = CompCor, 3 = Compcor and Global mean)');
        end

        if ~isfield(switches, 'skipfirst');
            switches.skipfirst = 0;
        end %skip this many frames at beginning of each run

        if ~isfield(switches, 'FDcut');
            switches.FDcut = input('Remove frames above FD (recommend 0.04, set to 100 to ignore)');
        end % threshold from petersen values to cut (recommend 0.04)

        if ~isfield(switches, 'DVARcut');
            switches.DVARcut = 1.05;
        end % percent above median to cutoff. zero means no cutoff (recommend 1.05)

        if ~isfield(switches, 'permutemask');
            switches.permutemask = 0;
        end % randomly scrambles the tmask

        if ~isfield(switches, 'figs');
            switches.figs = input('Create a tons of figs to QC your processing? (0=no, 1=yes)');
        end % generate figs to QC processing

        if ~isfield(switches, 'minframes');
            switches.minframes = input('Minimum frames for run exclusion? (recommend 400)');
        end % (recommend 400, unless you are afraid to throw away 1/4 subjects)

        if ~isfield(switches, 'TR');
            switches.TR = input('TR? (0.720 for HCP rfMRI data)');
        end % Acquisition rate in seconds

        if ~isfield(switches, 'FIX');
            switches.FIX = input('Use data run through FIX denoisifier? (0=no, 1=yes)');
        end

        if ~isfield(switches, 'corrtype');
            switches.corrtype = 'pearson';
        end % 'pearson' or 'partial'

        if ~isfield(switches, 'restart');
            switches.restart = input('Restart processing? (0=no (pick up where we left off), 1=yes)');
        end


        %% name output

        numregs = 25 + 3 * switches.GSR;

        if switches.permutemask
            outfilename = sprintf(fullfile(work_dir, 'outputs', '%s_norm_motion%d_GSR%d_crubPERM_frames%d_AveRun%d_fix%d_%sFC_dtseries.mat'),...
                switches.DataType, numregs,switches.GSR,switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        elseif switches.runs
            outfilename = sprintf(fullfile(work_dir, 'outputs', '%s_norm_motion%d_GSR%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_runs_dtseries.mat'),...
                switches.DataType, numregs,switches.GSR,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        else
            outfilename = sprintf(fullfile(work_dir, 'outputs', '%s_norm_motion%d_GSR%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_dtseries.mat'),...
                switches.DataType, numregs,switches.GSR,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        end
        numregs = 25 + 3 * switches.GSR;



        %% restart

        if switches.restart == 0
            try
                load(outfilename)
            catch
                fprintf('Data has not been processed yet.')
            end
        end

        %% CAROLINA not doing FC for all grayordinates

        if switches.permutemask
            outfilename = sprintf(fullfile(work_dir, 'outputs', 'sub%s_%s_norm_motion%d_GSR%d_crubPERM_frames%d_AveRun%d_fix%d_%sFC_dtseries.mat'),...
                Subjlist{ll},switches.DataType, numregs,switches.GSR,switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        elseif switches.runs
            outfilename = sprintf(fullfile(work_dir, 'outputs', 'sub%s_%s_norm_motion%d_GSR%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_runs_dtseries.mat'),...
                Subjlist{ll},switches.DataType, numregs,switches.GSR,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        else
            outfilename = sprintf(fullfile(work_dir, 'outputs', '%s_norm_motion%d_GSR%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_dtseries.mat'),...
                switches.DataType, numregs,switches.GSR,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        end





        sub = num2str(Subjlist{ll});
        fprintf('Loading Subj: %s\n', sub)
        if ll ==  2
            fprintf('Here')
        end
        sub_dir = fullfile(source_dir, sub);

        try

            dt=[];
            if iTask==1
                RestReg=cell(1,NR);
            elseif iTask==2
                TaskReg=cell(1,numel(TaskCell));
            end
            for t = 1:NR % runs

                %if (iTask==1)&&switches.FIX % passed through ICA-FIX
                %    runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas_hp2000_clean']);
                %else
                runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas']);
                %end
                runname
                [runpath, runfile] = fileparts(runname);
                dtseries = fullfile([runname '.dtseries.nii']);
                dtseries
                % read image
                disp('Reading in image')
                raw_cii = ft_read_cifti_mod(dtseries); disp('Read');% this is Siegel's modified version. Wasn't provided at the beginning
                if iTask==1
                    dtmean = nanmean(raw_cii.data,2);
                    dtstd = nanstd(raw_cii.data, [], 2);
                    RestMean{t}=dtmean;
                    RestSTD{t}=dtstd;
                else
                    %% Just start with doing the same LR/RL using session 1
                    Task2Rest=t;
                    TaskMatch{iTask-1}{t}=Task2Rest;
                    TaskMean{ll}{iTask}{t} = nanmean(raw_cii.data,2);
                    TaskSTD{ll}{iTask}{t} = nanstd(raw_cii.data, [], 2);
                    dtmean=RestMean{Task2Rest};
                    dtstd=RestSTD{Task2Rest};
                end
                tp = size(raw_cii.data, 2);
                vox = size(raw_cii.data, 1);

                % DEMEAN
                raw_cii.data=raw_cii.data-repmat(dtmean,1,size(raw_cii.data,2));
                % VARIANCE NORMALIZE

                raw_cii.data = raw_cii.data./repmat(dtstd,1,size(raw_cii.data,2));
                %% Reintroduced Detrending: 6/24
                % DEDRIFT
                size(raw_cii.data)
                raw_cii.data = detrend(raw_cii.data')';

                % LOAD 24 PARAM NUISANCE REGRESSION
                MV = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors.txt'));
                MVdt = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors_dt.txt'));
                %% Don't zscore MV MVdt
                regs=[MV MVdt];
                %regs = zscore([MV MVdt]);
                X.regs(iTask).run(t).MV=MV;
                X.regs(iTask).run(t).MVdt=MVdt;


                if iTask~=1
                    X.EV{iTask}{t}=AddEVs(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'EVs'));
                end


                %% DONT NEED TO USE 24 REGS IN FIX-ICA
                if switches.GSR == 0

                    regs = [regs ones(size(regs,1),1)];

                elseif switches.GSR == 1 %(mean WM, CSF, GM timecourses)

                    GS = mean(raw_cii.data,2);
                    global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'WM_CSF_GM_regs.txt')));
                    global_regs_prime = [[0 0 0]; diff(global_regs)];
                    regs = [ global_regs global_regs_prime ones(size(regs,1),1)];
                    % ADD TEMPORAL DERIVATIVES OF REGS

                elseif switches.GSR == 2 % compcor

                    %                if switches.FIX
                    %                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'WM_CSF_CompCor_regs_hp200_clean.txt')));
                    %                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'/WM_CSF_CompCor_regs.txt')));
                    %                end

                    regs = [ CompCor_regs(:,1:10) ones(size(regs,1),1)];


                elseif switches.GSR == 3 % compcor and global mean?

                    %                if switches.FIX
                    %                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_CompCor_regs_hp200_clean.txt')));
                    %                    global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_GM_regs_hp200_clean.txt')));
                    %                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_CompCor_regs.txt')));
                    global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_GM_regs.txt')));
                    %               end

                    regs = [ CompCor_regs(:,1:10) global_regs(:,3) ones(size(regs,1),1)];
                end % regs

                if switches.GSR>=2
                    X.regs(iTask).run(t).CompCor=CompCor_regs(:,1:10);
                    if switches.GSR==3
                        X.regs(iTask).run(t).GSR=global_regs(:,3);
                    end
                end
                X.switches(ll).run(t)=switches;
                if iTask==1
                    X.RestReg{t}=regs;
                else
                    X.TaskReg{iTask-1}{t}=regs;
                end
                %% APPLY REGRESSORS
                %            disp('Applying Regressors')
                %            for r = 1:vox
                %                [B,~,R] = regress(raw_cii.data(:,r), regs);
                %                raw_cii.data(:,r) = R;
                %                betas(r,:) = B;
                %            end % regress
                %            disp('done with vox')
                QC(ll).run(t).tmask = zeros(size(raw_cii.time));

                %% FD Scrubbing

                if switches.FDcut

                    [status,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [ sub '_' tseries{t} '_FD.txt'])]);

                    QC(ll).run(t).FD=single(importdata(strcat(filename)));

                    QC(ll).run(t).tmask(1) = 1;
                    pfd = find(QC(ll).run(t).FD>switches.FDcut);
                    QC(ll).run(t).tmask(pfd) = 1;

                end

                %% DVAR Scrubbing
                disp('DVAR Scrubbing')


                [status,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [sub '_' tseries{t} '_DVARS.txt'])]);

                QC(ll).run(t).DV=single(importdata(strcat(filename)));

                subDVcut = median(QC(ll).run(t).DV(2:end))*switches.DVARcut ;
                QC(ll).run(t).DVcut=subDVcut;
                if switches.DVARcut
                    pdv = find(QC(ll).run(t).DV>subDVcut);
                    QC(ll).run(t).tmask(pdv)=1;
                    QC(ll).run(t).tmask=single(QC(ll).run(t).tmask(1:length(QC(ll).run(t).DV)));

                end % DVARS scrubbing

                tmask = QC(ll).run(t).tmask;

                %% permute
                if switches.permutemask %% OPTIONALLY PERMUTE TMASK

                    tmaskorder = randperm(length(QC(ll).run(t).tmask));
                    tmask = QC(ll).run(t).tmask(tmaskorder);

                end % permute tmask


                %% Apply mask
                disp('Applying Mask')
                dt = raw_cii.data;

                for iY=1:(nParcLevels)
                    Combo=ParcMask{iY};
                    ParcLabs=unique(Combo(Combo~=0)');
                    ParcDat=zeros(numel(ParcLabs),size(dt,2));
                    whos ParcLabs
                    whos ParcDat
                    whos dt
                    whos Combo
                    for iParc=1:numel(ParcLabs)
                        disp([


                        iParc iTask])
                        ParcDat(ParcLabs(iParc),:)=nanmean(dt(Combo==ParcLabs(iParc),:),1);
                    end
                    tt=tmask;
                    ParcDat(:,tt==1)=interp1(find(tt==0),ParcDat(:,tt==0)',find(tt==1))';
                    if strcmpi(Atlas(1),'G')
                        ParcDat(1:333,:)=reorganizeSRI(ParcDat(1:333,:),ParcMask{2});
                    end
                    FullDat{iY,iTask}{t}=ParcDat;
                end


                QC(ll).frames(t) = sum(tmask==0);
                QC(ll).excluded(1:NR) = 0;

            end %% series

            %% Some QC stuff?
            disp('QC-ing')
            if QC(ll).frames(1) < switches.minframes || QC(ll).frames(2) < switches.minframes
                %dt(1:2)=[];
                %QC(n).frames(1:2)=0;
                QC(ll).excluded(1:NR)=1;
            end

            excluded(ll) = 0; %MITCH Set default values to initiate variable

            if sum(QC(ll).excluded) == NR;
                excluded(ll)=1;
            end

        catch
            sprintf('Subject Failed.\n')
            bugs(ll) = 1;

        end

        if excluded(ll)==1
            sprintf('Subject Excluded.\n')
        end

        myQC=QC(ll);
        myQC.BaseName=outfilename;
        FullQC{iTask}=myQC;
        if iTask==1
            TaskName{iTask}='Rest';
        else
            TaskName{iTask}=TaskAll{iTask-1};
        end
    end
    X.Proc={'Av','Av','Av','Av','n','y','y'};
    X.GSR=GSR;
    X.QC=FullQC;
    X.Notes='Doesn`t use ICA-FIX or DVAR cuts to match task. Does use detrending';
    X.TaskName=TaskName;
    X.TaskMatch=TaskMatch;

    bDir='/scratch2/Singh/HCP_All/GSR';
    disp('Preparing to Save')
    if strcmpi(Atlas(1),'y')
        YeoSize={'01','02','04','06','08','10'};
        for iY=1:size(FullDat,1)
            ParcName=strcat('Y',YeoSize{iY});
            %disp(ParcName)
            X.Dat=FullDat(iY,:);
            save(strcat(bDir,num2str(GSR),'/sub',sub,ParcName,'_All','.mat'),'X');
        end
    else
        ParcName=Atlas(1);
        X.Dat=FullDat;
        save(strcat(bDir,num2str(GSR),'/sub',sub,ParcName,'_All','.mat'),'X');
    end
end
end