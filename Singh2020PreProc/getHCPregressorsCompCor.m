%% Modified by Ruiqi Chen
%
%% getHCPregressoresCompCor_final
% The purpose of this script is to use freesurfer segmentations as well as
% volume rfMRI data to generate both average tissue regressors as well as
% multiple component-based tissue regressors (as in Behzadi, 2007).
%
% (c) 2016 Washington University in Saint Louis
% Author: Joshua S. Siegel, jssiegel@wustl.edu
% All Rights Reserved
function getHCPregressorsCompCor(Subject, in_dir, out_dir, tseries, switches)

% FIX or non-FIX
if switches.FIX
    suffix = '_hp2000_clean';
else
    suffix = '';
end

% Handle missing data
missFlag = zeros(numel(tseries), 1);
for i = 1:numel(tseries)
    if ~isfile(fullfile(in_dir, Subject, 'MNINonLinear', 'Results', tseries{i}, [tseries{i} suffix '.nii.gz']))
        missFlag(i) = 1;
    end
end
missFlag = logical(missFlag);
tseries = tseries(~missFlag);
if isempty(tseries)
    return
end

%%GMkeeprgns	= [ 7 8 10 11 12 13 16 17 18 26 28 46 47 49 50 51 52 53 54 58 60 ];

compcorregs = 1;
savemeanregs = 1;
figs = 1;

%%matlabpool(3)
%Make structure elements for white matter and CSF
seWM = strel('disk',2,0)
seCSF = strel('disk',1,0)

WMkeeprgns	= [ 2 41 ];%Keep the Areas in the aparc+aseg.nii.gz that have one of these values
CSFkeeprgns	= [ 4 14 15 43 ];%are these not present
tic

sub= Subject
asegname = fullfile(in_dir, sub, '/MNINonLinear/aparc+aseg.nii.gz'); 
if exist(asegname,'file')
    try
        aseg = load_nifti(asegname);

        % make WM/GM/CSF masks
        WMmask = ismember(aseg.vol,WMkeeprgns);
        CSFmask = ismember(aseg.vol,CSFkeeprgns);

        rfMRIname = fullfile(in_dir, sub, '/MNINonLinear/Results', tseries{1}, [tseries{1} suffix '.nii.gz']);
        TC = load_nifti(rfMRIname);%I load in the Tseries here to pull the dimensions for resampling

        mask = [];
        % downsample masks
        %MITCH: I modified these next lines to account for DMCC size images
        WMmask = resample3Dimage(single(WMmask),aseg.dim(2)/TC.dim(2));%resamle the mask to that it will line up with the 
        WM2=imerode(WMmask==1,seWM);%eroded white matter mask
        mask(:,1) = WM2(:);
        CSFmask = resample3Dimage(single(CSFmask),aseg.dim(2)/TC.dim(2));%resample csf mask
        %CSF2=imerode(CSFmask==1,seCSF);
        mask(:,2) = CSFmask(:)==1;%store Csf mask as vector 
        mask = logical(mask);%convert 1's & 0's to logical values 
        if figs
            figure; subplot(2,2,1);imagesc(WMmask(:,:,46)==1);title('original image');hold on;
            subplot(2,2,3);imagesc(WM2(:,:,46));title('eroded image');hold on;
            subplot(2,2,2);imagesc(CSFmask(:,:,48)==1);title('original image');hold on;
        %      subplot(2,2,4);imagesc(CSF2(:,:,48));title('eroded image');
        end

        if savemeanregs
            GMmask = (aseg.vol>1000);%only keep values from aseg that are greater than 1000
                %MITCH: I modified these next lines to account for DMCC size images 
            GMmask = resample3Dimage(single(GMmask),aseg.dim(2)/TC.dim(2));%resample the image
            mask(:,3) = GMmask(:)==1;%Store as a vector in the mask
            %mask(mask<1)=0;
        end

        for t=1:length(tseries)

            outputname = fullfile(out_dir, sub, 'Results', tseries{t}, ['WM_CSF_CompCor_regs' suffix '.txt']);
            
            if ~exist(outputname,'file')
                
                rfMRIname = fullfile(in_dir, sub, '/MNINonLinear/Results', tseries{t}, [tseries{t} suffix '.nii.gz']);
                TC = load_nifti(rfMRIname);
                
                disp(size(TC.vol))
                    %MITCH: I modified these next lines to account for DMCC size images
                %save this dimension for later
                dimension = TC.dim(5)
                TC=single(reshape(TC.vol,TC.dim(2)*TC.dim(3)*TC.dim(4),TC.dim(5))); %Change into a vector
                % apply masks
                if savemeanregs
                    disp('Mask Size and TC size')
                    size(mask)
                    size(TC)
                    meanregs = TC'*mask;
                    meanregs = zscore(meanregs);
                end
                % PCA
                regs = [];
                for m=1:2
                    %% Demean, detrend, variance normalize
                    R_TC=TC(mask(:,m),:);
                    R_TC = detrend(R_TC')';
                    ptstd = nanstd(R_TC,[],2);
                    %R_TC = R_TC./repmat(ptstd,1,size(R_TC,2));
                    
                    %% PCA
                    [coeff, score, latent] = pca(R_TC);
                    regs = [regs coeff(:,1:5)];
                end
                %pve=cumsum(latent)./sum(latent);
                if compcorregs
                    fid = fopen(outputname,'w');
                    for i=1:size(regs,1)-1
                        fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',regs(i,:));
                    end
                    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f',regs(dimension,:));
                    fclose(fid);
                end
                
                if savemeanregs
                    fid = fopen(fullfile(out_dir, sub, 'Results', tseries{t}, ['WM_CSF_GM_regs' suffix '.txt']),'w');
                    for i=1:size(meanregs,1)-1
                        fprintf(fid,'%.4f\t%.4f\t%.4f\n',meanregs(i,:));
                    end
                    fprintf(fid,'%.4f\t%.4f\t%.4f',meanregs(dimension,:));
                    fclose(fid);
                end
            end
            toc
        end
        fprintf('Successfully generated regressors for subject %s.\n', sub)
    catch
        fprintf('Unknown error when generating regressors for subject %s.\n', sub)
    end
end
toc
%matlabpool close