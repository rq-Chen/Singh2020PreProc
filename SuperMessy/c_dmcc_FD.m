function c_dmcc_FD(in_dir,out_dir)

% For DMCC data in nil-bluearc 

% example:
% in_dir='/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/MINIMALLY_PREPROCESSED';
% out_dir = '/data/nil-bluearc/ccp-hcp/scratch/motion_timeseries/temp_FD'

subjects = strsplit(ls(in_dir));
for i=1:length(subjects)

    if ~strcmp(subjects{i}, '') && ~strcmp(subjects{i}, '-i')
        aux_dir     = fullfile(in_dir, subjects{i}, 'MNINonLinear', 'Results');
        bolds  = strsplit(ls(aux_dir));
        
        for j = 1:length(bolds)
            if ~isempty(strfind(bolds{j}, 'rfMRI_Rest')) || ~isempty(strfind(bolds{j}, 'tfMRI'))
                movregs = fullfile(aux_dir, bolds{j}, 'Movement_Regressors.txt');
                get_FD(subjects{i}, bolds{j}, movregs, out_dir)
            end
        end
        
    end 
end
end
