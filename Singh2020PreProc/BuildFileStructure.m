%build the File Structure to store the result data.
%This should make a results folder in the out_dir for each subject and a folder
%for each tseries. this should also make an 'outputs' folder in the out_dir
%where the .mat data will be stored
%may need to be edited if we decide to store the results with the preproccessed data this. 

function BuildFileStructure(workDir, Subject, tseries)
    if ~exist(fullfile(workDir, Subject), 'dir');
        mkdir(fullfile(workDir, Subject))
    end
    if ~exist(fullfile(workDir, Subject,'Results'), 'dir');
        mkdir(fullfile(workDir, Subject, 'Results'))
    end
    if ~exist(fullfile(workDir, 'outputs'), 'dir');
        mkdir(fullfile(workDir, 'outputs'))
    end 
    for n=1:length(tseries)
        if ~exist(fullfile(workDir, Subject,'Results', tseries{n}),'dir');
            mkdir(fullfile(workDir, Subject,'Results', tseries{n}));
            disp(fullfile(workDir, Subject, 'Results', tseries{n}))
        end
    end
end
