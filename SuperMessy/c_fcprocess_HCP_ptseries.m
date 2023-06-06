function c_fcprocess_HCP_ptseries(switches, Subjlist, source_dir, work_dir, tseries, Parc)

%TODO[c]: fix the paths mess --> as input
% Make GL 32k parcellated connectivity


%% initialize variables

       
NR = length(tseries);   % number of runs (LR and RL in HCP)

%% Load subject list

%%load(fullfile(work_dir, 'Subjlist.mat')) % TODO: as input too?
%% Load a bunch of settings (You can preset them in the switches struct)
%This Just makes sure that all the switches are set
if ~isfield(switches, 'AveRun');
    switches.AveRun = input('Do you want to concatenate timecourses (0) or calc FC by run and then average (1)?');
end

if ~isfield(switches, 'GSR');
    switches.GSR = input('Timecourse Regression? (0 = No, 1 = GSR, 2 = CompCor)');
end

if ~isfield(switches, 'bpss');
    switches.bpss = input('Bandpass filter? (0=no, 1=yes)');
end

if ~isfield(switches, 'bpss_lo') && switches.bpss == 1;
    switches.bpss_lo = input('Low pass? (recommend 0.08)');
end % (recommend 0.08)

if ~isfield(switches, 'bpss_hi') && switches.bpss == 1;
    switches.bpss_hi = input('Hi pass? (recommend 0.009)');
end % (recommend 0.009)

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

%% Set up bandpass butterworth filter (CAROLINA: normalized frequency!!!)
%Mitch: I changes this to account for alternative TRs
filterorder = 1;
lopass = switches.bpss_lo / (0.5/switches.TR);
hipass = switches.bpss_hi / (0.5/switches.TR);
[butta, buttb] = butter(filterorder, [hipass lopass]);



%% name output
%Names the .mat output that contains all the important info
numregs = 25 + 3 * switches.GSR;

if switches.permutemask
    outfilename = sprintf(fullfile(work_dir,'outputs', '%s__%s_motion%d_GSR%d_bpss%d_scrubPERM_frames%d_AveRun%d_fix%d_%sFC.mat'),...
        switches.DataType,Parc.name,numregs,switches.GSR,switches.bpss,switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);

elseif switches.runs
    outfilename = sprintf(fullfile(work_dir, 'outputs', '%s__%s_norm_motion%d_GSR%d_bpss%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC_runs.mat'),...
        switches.DataType,Parc.name,numregs,switches.GSR,switches.bpss,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);

else
    outfilename = sprintf(fullfile(work_dir,'outputs', '%s__%s_norm_motion%d_GSR%d_bpss%d_scrub%d_frames%d_AveRun%d_fix%d_%sFC.mat'),...
        switches.DataType,Parc.name,numregs,switches.GSR,switches.bpss,100*(switches.FDcut),switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);

end


%% restart

if switches.restart == 0
    try
        load(outfilename)
    catch
        fprintf('Data has not been processed yet.')
    end
end

if ~exist('HCP_GL_FC','var')
    HCP_GL_FC = nan(length(Subjlist), Parc.NP * (Parc.NP-1)/2); 
    % HCP_GL_FC = nan(length(Subjlist),324*323/2); % CAROLINA: THESE VALUES
    % WERE FOR THE 324 GORDON PARCELS    
end

if switches.runs
%   run_FC = nan(length(Subjlist),4,324*323/2);
    run_FC = nan(length(Subjlist), NR, Parc.NP * (Parc.NP-1)/2);
end


excluded = zeros(length(Subjlist),1);
bugs = zeros(length(Subjlist),1);
dolist = find(isnan(sum(HCP_GL_FC,2)));
Subjlist = Subjlist(dolist);


QC = [];
regs = [];
runs = zeros(length(Subjlist),1);
switches


%% Do the stuff

for ll = 1:length(Subjlist) % subjects

    n = dolist(ll);
    sub = num2str(Subjlist{ll});
    fprintf('Loading Subj: %s\n', sub)
    
    sub_dir = fullfile(source_dir, sub);

    try
        
        pt=[];

        for t = 1:NR % runs 
            
            if switches.FIX %passed through ICA-FIX
                runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas_hp2000_clean']);    
            else
                runname = fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t}, [tseries{t} '_Atlas']);           
            end

            
            %% Parcellation TODO[carolina]: it would be better if the parcellation was already done outside this
                         
            [runpath, runfile] = fileparts(runname);
            ptseries = fullfile(work_dir, sub, 'Results', tseries{t}, [runfile '_' Parc.name '.ptseries.nii']);
            
            if ~exist(ptseries,'file')  % parcellate
                
                tic
                system(['wb_command -cifti-parcellate ' runname '.dtseries.nii ' Parc.dir ' COLUMN ' ptseries]);
                toc
                
            end
            
            % read image
            disp('Reading in image')
            
            raw_cii = ft_read_cifti_mod(ptseries); % this is Siegel's modified version. Wasn't provided at the beginning
            
            if ~Parc.ordered%if the Parcellation is not ordered then reorder it 
               raw_cii.data = reorganize(raw_cii.data, Parc);
            end
            
            ptmean = nanmean(raw_cii.data,2); %find the mean value of all the parcels across time 
            ptstd = nanstd(raw_cii.data, [], 2);%find the std of all parcels across time
            tp = size(raw_cii.data, 2); %find the number of time points 
            vox = size(raw_cii.data, 1);%find the number of voxels
            
            % DEMEAN
            % VARIANCE NORMALIZE
            raw_cii.data = raw_cii.data./repmat(ptstd,1,size(raw_cii.data,2));
            % DEDRIFT
            raw_cii.data = detrend(raw_cii.data');
            rawTC = raw_cii.data';
            
            % LOAD 24 PARAM NUISANCE REGRESSION
            MV = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors.txt'));
            MVdt = importdata(fullfile(sub_dir, 'MNINonLinear', 'Results', tseries{t},'Movement_Regressors_dt.txt'));
            regs = zscore([MV MVdt]);
            
            %% DONT NEED TO USE 24 REGS IN FIX-ICA PIPELINE
            if switches.GSR == 0
                
                regs = [regs ones(size(regs,1),1)];
                
            elseif switches.GSR == 1 %(mean WM, CSF, GM timecourses) file WM_CSF_GM_regs.txt needs to exist already!!
                
                GS = mean(raw_cii.data,2);
                global_regs = zscore(importdata(fullfile(work_dir, sub, 'Results', tseries{t},'WM_CSF_GM_regs.txt')));
                global_regs_prime = [[0 0 0]; diff(global_regs)];
                regs = [ global_regs global_regs_prime ones(size(regs,1),1)];
                % ADD TEMPORAL DERIVATIVES OF REGS
                
            elseif switches.GSR == 2
                
                if switches.FIX
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t},'WM_CSF_CompCor_regs_hp200_clean.txt')));
                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t},'/WM_CSF_CompCor_regs.txt')));
                end
                
                regs = [ CompCor_regs(:,1:10) ones(size(regs,1),1)];
 
                
            elseif switches.GSR == 3 %% (CompCor)
                
                if switches.FIX
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t}, 'WM_CSF_CompCor_regs_hp200_clean.txt')));
                    global_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t}, 'WM_CSF_GM_regs_hp200_clean.txt')));
                else
                    CompCor_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t}, 'WM_CSF_CompCor_regs.txt')));
                    global_regs = zscore(importdata(fullfile(work_dir, sub,'Results', tseries{t}, 'WM_CSF_GM_regs.txt')));
                end
                
                regs = [ CompCor_regs(:,1:10) global_regs(:,3) ones(size(regs,1),1)];
            end % regs
            
            %% APPLY REGRESSORS
            disp('Applying Regressors')
            for r = 1:vox
                
                [B,~,R] = regress(raw_cii.data(:,r), regs);
                raw_cii.data(:,r) = R;
                betas(r,:) = B;
                
            end % regress
            
            QC(n).run(t).tmask = zeros(size(raw_cii.time));
            
            %% FD Scrubbing
            disp('FD Scrubbing')
            if switches.FDcut                
                
                if exist('FDfilt','var')
                    
                    QC(n).run(t).FD=single(FDfilt(n,t,:));                    
                
                else
                    
                    % TODO[carolina] : create the filtered timeseries
                    %[status,filename] = system(['ls /data/cn4/collaborators/HCP/fcprocessing/fd_filt_266/*_' sub '_' tseries{t} '_FD_FILT.txt']);
                    [status,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [ sub '_' tseries{t} '_FD_FILT.txt'])]);
                    
                    % CAROLINA: This should be status ~= 0
                    % if status==1 (original version)
%                     if status ~=0
%                         % CAROLINA: what's this path here? I think is an
%                         % alternative search path for some of their data??
%                         % We don;t ahve this, so we won't care....
%                         [status,filename] = system(['ls /data/cn4/collaborators/HCP/fcprocessing/FD/*_' sub '_' tseries{t} '_FILT.txt']);
%                         notfound(ll,t)=1;
%                     end
                    
                    QC(n).run(t).FD=single(importdata(strcat(filename)));
                    
                end % FD scrubbing
                
                QC(n).run(t).tmask(1) = 1;
                pfd = find(QC(n).run(t).FD > switches.FDcut); % indexes of FD censored frames
                QC(n).run(t).tmask(pfd) = 1;
                
            end
            
            %% DVAR Scrubbing
            disp('DVAR Scrubbing')
            if switches.DVARcut
                % Calc DV on parcels
                %                 dcii=diff(raw_cii.data);
                %                 DVparcel=[0 ; sqrt(mean(dcii.^2,2))];
                
                [status,filename] = system(['ls ' fullfile(work_dir, sub, 'Results', tseries{t}, [sub '_' tseries{t} '_DVARS_FILT.txt'])]);
                
                % CAROLINA: This should be status ~= 0
                % if status==1 (original version)
                % same as for FDs
                if status ~=0
                    [status,filename] = system(['ls /data/cn4/collaborators/HCP/fcprocessing/DVARS/*_' sub '_' tseries{t} '_DV_FILT.txt']);
                end
                
                QC(n).run(t).DV=single(importdata(strcat(filename)));
                
                subDVcut = median(QC(n).run(t).DV(2:end))*switches.DVARcut ;
                pdv = find(QC(n).run(t).DV > subDVcut);     % indexes of DVARS censored frames
                QC(n).run(t).tmask(pdv)=1;
                QC(n).run(t).tmask=single(QC(n).run(t).tmask(1:length(QC(n).run(t).DV)));
                
            end % DVARS scrubbing
            
            tmask = QC(n).run(t).tmask;
            
            %% permute
            if switches.permutemask %% OPTIONALLY PERMUTE TMASK
                
                tmaskorder = randperm(length(QC(n).run(t).tmask));
                tmask = QC(n).run(t).tmask(tmaskorder);
                
            end % permute tmask
            
            %% INTERPOLATE MISSING FRAMES
            if switches.bpss && switches.FDcut
                
                cuts = find(tmask==1);
                k = 1;
                C{k} = cuts(k);
                
                for j = 2:length(cuts)
                    td = cuts(j)-cuts(j-1);
                    if td == 1
                        C{k} = [C{k} cuts(j)];
                    else
                        k = k + 1;
                        C{k} = cuts(j);
                    end
                end
                
                for m = 1:length(C)
                    
                    k = C{m}(end) - C{m}(1)+1;
                    fill = nan(k,vox);
                    if m == 1
                        before = zeros(1,vox);
                    else
                        before = raw_cii.data(C{m}(1)-1,:);
                    end
                    if C{m}(end) == tp
                        after = zeros(1,vox);
                    else
                        after = raw_cii.data(C{m}(end)+1,:);
                    end
                    for j = 1:k
                        fill(j,:) = (j*after+k*before)/(j+k);
                    end
                    raw_cii.data(C{m}(1):C{m}(end),:) = fill;
                    
                end
                
            end % interpolate missing frames
            
            %% BANDPASS
            disp('Running Bandpass')
            if switches.bpss
                
                prefilt = raw_cii.data;
                pad = ceil(1/switches.bpss_hi);
                prefilt = [zeros(pad,vox) ; prefilt ; zeros(pad,vox)];
                
                % apply butterworth bpass filter
                postfilt = filtfilt(butta,buttb,prefilt);
                prefilt = prefilt((pad+1):(end-pad),:);
                postfilt = postfilt((pad+1):(end-pad),:);
                raw_cii.data = postfilt;
                
            end % bandpass
                       
            
            %% Apply mask
            disp('Apply mask')
            pt{t} = raw_cii.data';
            pt{t} = pt{t}(:,tmask==0);
            
            QC(n).frames(t) = sum(tmask==0);
            QC(n).excluded(1:NR) = 0;
            
            %% FIGS
            if switches.figs % Make some QC figs
                
                % check spectrum
                for i=1:vox
                    [P,F]=pwelch(raw_cii.data(:,i),50,[],256,1/switches.TR);
                    P_runs(i,:)=P;
                end
                
                meanP=nanmean(P_runs);
                steP=std(P_runs)/sqrt(size(P_runs,1));
                figure
                boundedline(F,meanP,steP,'cmap',[0,0,0]); hold on
                xlim([0.001 0.2]);
                xlabel('Frequency');ylabel('Power')
                
                
                % check timecourse
                NUS = [zscore(QC(n).run(t).DV) zscore(QC(n).run(t).FD)];
                figure;
                ax1 = subplot(5,1,1)
                plot(1:length(raw_cii.time),QC(n).run(t).FD,'r')
                hold on;line([0 length(raw_cii.time)],[switches.FDcut switches.FDcut])
                ax2 = subplot(5,1,2)
                plot(1:length(raw_cii.time),QC(n).run(t).DV)
                hold on;line([0 length(raw_cii.time)],[subDVcut subDVcut])
                ylim([nanmean(QC(n).run(t).DV)*.8 nanmean(QC(n).run(t).DV)*1.2])
                ax3 = subplot(5,1,3)
                imagesc(rawTC)
                colormap('gray');
                title('Detrended only')
                linkaxes([ax1,ax2],'x')
                ax4 = subplot(5,1,4)
                imagesc(postfilt')
                colormap('gray')
                title('filtered')
                raw_cii.data(tmask==1,:) = nan;
                ax5 = subplot(5,1,5)
                imagesc(raw_cii.data')
                colormap('gray')
                title('final')
                linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
                
                % check FC
                figure
                GLvisFC(atanh(corr(pt{t}')),[-0.4 1.2],['run ' num2str(t)]); % CAROLINA: we don't have this function
                
                
            end % figures
            
        end %% series
        %% Some QC stuff?

        % CAROLINA : comented some QC stuff, because we don't have other
        % files...( we just have 2 runs per experiment)
        
%         QC(n).frames(1) < switches.minframes || QC(n).frames(2) < switches.minframes % hard coded for number of runs 
%             %pt(1:2)=[];
%             %QC(n).frames(1:2)=0;
%             QC(n).excluded(1:NR)=1;
%         end
        %MITCH: expanded to allow for variable size run
        QC(n).excluded = QC(n).frames < switches.minframes
        
        %% CORRELATIONS
        disp('Running Correlations')
        if switches.AveRun % PEARSON CORRELATION ON TIMESERIES BY RUN
            
            for i = 1:NR;
                if QC(n).excluded(i) == 0
                    switch switches.corrtype
                        case 'pearson'
                            RunCorr(i,:) = ExtractDataAboveDiagonal(atanh(corr(pt{i}')));
                        case 'partial'
                            RunCorr(i,:) = ExtractDataAboveDiagonal(nets_netmats(pt{i}',1,'ridgep',0.01));
                    end
                end
            end
            
            if switches.runs
                QC(n).FC = RunCorr; % CAROLINA: error if all runs were excluded, since RunCorr doesn't exist
            end
            
            SeedCorrMatrix = squeeze(mean(RunCorr(QC(n).excluded==0,:),1));
        
        else % CONCATENATE WHOLE TIMESERIES THEN CORRELATE

            TC = [];
            for i = 1:NR
                if QC(n).excluded(i) == 0
                    TC = [TC pt{i}];
                end
            end
            
            switch switches.corrtype
                case 'pearson'
                    disp('Correlating')
                    SeedCorrMatrix = ExtractDataAboveDiagonal(atanh(corr(TC')));
                case 'partial'
                    SeedCorrMatrix = ExtractDataAboveDiagonal(nets_netmats(TC',1,'ridgep',0.01));
            end
        end % average or concatenate results
        
        SingleSubOut = sprintf(fullfile(work_dir,sub, '%s__%s__%s_motion%d_GSR%d_bpss%d_scrubPERM_frames%d_AveRun%d_fix%d_%sFC.mat'),...
        sub,switches.DataType,Parc.name,numregs,switches.GSR,switches.bpss,switches.minframes,switches.AveRun,switches.FIX,switches.corrtype);
        subQC = QC(n)
        subBugs = bugs(n)
        save(SingleSubOut,'SeedCorrMatrix','excluded','subBugs','subQC','switches','Parc','-v7.3')

        HCP_GL_FC(n,:) = SeedCorrMatrix;
        
        if sum(QC(n).excluded) == NR;
            excluded(n)=1;
        end
        
    catch
        sprintf('Subject Failed.\n')
        bugs(n) = 1;
        
    end
    
    if excluded(ll)==1
        sprintf('Subject Excluded.\n')
    end
        
    
end

%save the data for all the subjects ran
save(outfilename,'HCP_GL_FC','excluded','bugs','QC','switches','Parc','Subjlist','-v7.3')

%this wil rebuild the matrix to allow for visualizaion of the unwrapped
%cross-run averaged data
for ll = 1:length(Subjlist)
    figure('Name', ['Cross Run Averaged FC of ' switches.DataType 's ' Subjlist{ll} ])
    n = Parc.NP;
    B = zeros(n);
    [ii jj] = ndgrid(1:n);
    B(ii>jj) = HCP_GL_FC(ll, :).';
    B = B.' + B;
    imagesc(B);
    colorbar;
end

end



