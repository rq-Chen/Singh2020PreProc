% For hcp data: Movement Reegressors are in degrees, not radians

function c_calc_FDretest(in_dir, out_dir, subjects, tseries)
    
    for i=1:length(subjects)
        %if the strings are not empty or -i
        if ~strcmp(subjects{i}, '') && ~strcmp(subjects{i}, '-i')
            %put the file together <indir>/<subject>/MNINonlinear/Results
           
            
            %for the files in results
            for j = 1:length(tseries)
                 aux_dir = fullfile(in_dir, strcat(subjects{i},'_3T/RESOURCES/',tseries{j},'_FIX'),tseries{j});% 'MNINonLinear', 'Results');
                 disp(tseries{j})
                %if there are LR results or RL results
                if ~isempty(strfind(tseries{j}, 'AP')) || ~isempty(strfind(tseries{j}, 'PA')) ...
                   || ~isempty(strfind(tseries{j}, 'LR')) || ~isempty(strfind(tseries{j}, 'RL'))
                    %store the full file name in movregs
                    movregs = fullfile(aux_dir, tseries{j}, 'Movement_Regressors.txt');
                    disp(movregs)
                    %if movreg is a file
                    if exist(movregs, 'file')
                        %get FD
                        get_FD(subjects{i}, tseries{j}, movregs, out_dir)
                    end
                end
            end
        
        end
    end
end


