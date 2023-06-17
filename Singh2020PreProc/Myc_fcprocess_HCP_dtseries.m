function Myc_fcprocess_HCP_dtseries(switches, Subjlist, source_dir, work_dir, tseries)
%TODO[c]: fix the paths mess --> as input
% Make GL 32k parcellated connectivity


%% initialize variables
%TODO : series names as input too?
NR = length(tseries);   % number of runs (LR and RL in HCP)


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

   

% %% restart

% if switches.restart == 0
%     try
%         load(outfilename)
%     catch
%         fprintf('Data has not been processed yet.')
%     end
% end

%% CAROLINA not doing FC for all grayordinates

switches

%% Do the stuff
for ll = 1:length(Subjlist) % subjects
    
    if switches.permutemask
        outfilename = sprintf(fullfile(work_dir, 'outputs', ...
            'sub%s_%s_norm_motion%d_GSR%d_crubPERM_frames%d_AveRun%d_fix%d_%sFC_dtseries.mat'), ...
            Subjlist{ll},switches.DataType, numregs,switches.GSR,switches.minframes,...
            switches.AveRun,switches.FIX,switches.corrtype);
    elseif switches.runs
        outfilename = sprintf(fullfile(work_dir, 'outputs', ...
            'sub%s_%s_norm_motion%d_GSR%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_runs_dtseries.mat'),...
            Subjlist{ll},switches.DataType, numregs,switches.GSR,100*(switches.FDcut),switches.minframes,...
            switches.AveRun,switches.FIX,switches.corrtype);
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
            
    dt = cell(1, NR);
    myQC = struct();
    myQC.run = struct('regs', cell(1, NR), 'tmask', cell(1, NR), 'FD', cell(1, NR), 'DV', cell(1, NR));
    myQC.frames = nan(1, NR);  

    for t = 1:NR % runs 

        try    
            if switches.FIX % passed through ICA-FIX
                runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas_MSMAll_hp2000_clean']);
            else
                runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas']);
            end

            dtseries = fullfile([runname '.dtseries.nii']);
            
            % read image
            disp('Reading in image')
            raw_cii = ft_read_cifti_mod(dtseries); % this is Siegel's modified version. Wasn't provided at the beginning
            % dtmean = nanmean(raw_cii.data,2);
            dtstd = nanstd(raw_cii.data, [], 2);
            % tp = size(raw_cii.data, 2);
            vox = size(raw_cii.data, 1);
            
            % DEMEAN
            % VARIANCE NORMALIZE
            raw_cii.data = raw_cii.data./repmat(dtstd,1,size(raw_cii.data,2));
            % DEDRIFT
            raw_cii.data = detrend(raw_cii.data');
            % rawTC = raw_cii.data';
            
            % LOAD 24 PARAM NUISANCE REGRESSION
            MV = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors.txt'));
            MVdt = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors_dt.txt'));
            regs = zscore([MV MVdt]);
            
            %% DONT NEED TO USE 24 REGS IN FIX-ICA            
            if switches.GSR == 0
                
                regs = [regs ones(size(regs,1),1)];
                
            elseif switches.GSR == 1 %(mean WM, CSF, GM timecourses)
                
                % GS = mean(raw_cii.data,2);
                global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'WM_CSF_GM_regs.txt')));
                global_regs_prime = [[0 0 0]; diff(global_regs)];
                regs = [ global_regs global_regs_prime ones(size(regs,1),1)];
                % ADD TEMPORAL DERIVATIVES OF REGS

            elseif switches.GSR == 2 % compcor
                
                if switches.FIX
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'WM_CSF_CompCor_regs_hp200_clean.txt')));
                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'/WM_CSF_CompCor_regs.txt')));
                end
                
                regs = [ CompCor_regs(:,1:10) ones(size(regs,1),1)];

                
            elseif switches.GSR == 3 % compcor and global mean?
                
                if switches.FIX
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_CompCor_regs_hp200_clean.txt')));
                    global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_GM_regs_hp200_clean.txt')));
                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_CompCor_regs.txt')));
                    global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t}, 'WM_CSF_GM_regs.txt')));
                end
                
                regs = [ CompCor_regs(:,1:10) global_regs(:,3) ones(size(regs,1),1)];
            end % regs
            
            myQC.run(t).regs=regs;
            
            %% APPLY REGRESSORS
            disp('Applying Regressors')
            for r = 1:vox
                
                [B,~,R] = regress(raw_cii.data(:,r), regs);
                raw_cii.data(:,r) = R;
                betas(r,:) = B;
                
            end % regress
            disp('done with vox')
            myQC.run(t).tmask = zeros(size(raw_cii.time));
            
            %% FD Scrubbing
            
            if switches.FDcut                
                
                [~,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [ sub '_' tseries{t} '_FD.txt'])]);
                                        
                myQC.run(t).FD=single(importdata(strcat(filename)));
                            
                myQC.run(t).tmask(1) = 1;
                pfd = find(myQC.run(t).FD>switches.FDcut);
                myQC.run(t).tmask(pfd) = 1;
                
            end
            
            %% DVAR Scrubbing
            disp('DVAR Scrubbing')
            if switches.DVARcut
                            
                [~,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [sub '_' tseries{t} '_DVARS.txt'])]);
                                            
                myQC.run(t).DV=single(importdata(strcat(filename)));
                
                subDVcut = median(myQC.run(t).DV(2:end))*switches.DVARcut ;
                pdv = find(myQC.run(t).DV>subDVcut);
                myQC.run(t).tmask(pdv)=1;
                myQC.run(t).tmask=single(myQC.run(t).tmask(1:length(myQC.run(t).DV)));
                
            end % DVARS scrubbing
            
            tmask = myQC.run(t).tmask;
            
            %% permute
            if switches.permutemask %% OPTIONALLY PERMUTE TMASK
                
                tmaskorder = randperm(length(myQC.run(t).tmask));
                tmask = myQC.run(t).tmask(tmaskorder);
                
            end % permute tmask
        
                        
            %% Apply mask
            disp('Applying Mask')
            dt{t} = raw_cii.data';
            dt{t} = dt{t}(:,tmask==0);
            
            myQC.frames(t) = sum(tmask==0);
    
        catch
            warning(['Error in subject ' sub ' session ' tseries{t}])
        end        
                    
    end

    if all(cellfun(@isempty, dt))
        warning(['All sessions failed for subject ' sub '! Will not save dtseries.'])
        continue
    end

    % %% Some QC stuff?
    % disp('QC-ing')
    % myQC.excluded(1:NR) = 0;
    % if myQC.frames(1) < switches.minframes || myQC.frames(2) < switches.minframes
    %     %dt(1:2)=[];
    %     %QC(n).frames(1:2)=0;
    %     myQC.excluded(1:NR)=1;
    % end
    
    % excluded = 0; %MITCH Set default values to initiate variable
    
    % if sum(myQC.excluded) == NR
    %     excluded=1;
    % end
    
    % if excluded==1
    %     sprintf('Subject Excluded.\n')
    % end

    save(outfilename,'dt', 'myQC', 'switches','-v7.3')

        
end
end
