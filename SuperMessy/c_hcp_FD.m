% For hcp data: Movement Reegressors are in degrees, not radians

function c_calc_FD(in_dir, out_dir)
    subjects = strsplit(ls(in_dir));

    for i=1:length(subjects)
    
        if ~strcmp(subjects{i}, '') && ~strcmp(subjects{i}, '-i')
            aux_dir     = fullfile(in_dir, subjects{i}, 'MNINonLinear', 'Results');
            bolds  = strsplit(ls(aux_dir));
        
            for j = 1:length(bolds)
                if ~isempty(strfind(bolds{j}, 'LR')) || ~isempty(strfind(bolds{j}, 'RL'))
                    movregs = fullfile(aux_dir, bolds{j}, 'Movement_Regressors.txt');
                    if exist(movregs, 'file')
                        get_FD(subjects{i}, bolds{j}, movregs, out_dir)
                    end
                end
            end
        
        end
    end
end


