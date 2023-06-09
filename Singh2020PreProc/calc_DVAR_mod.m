function calc_DVAR_mod(in_dir, out_dir, subjlist, tseries)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modified by Ruiqi Chen
    % Whether to use DVARS.sh or the MATLAB package "DVARS". See README.
    use_old_fslDVARS = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    script_dir = fileparts(mfilename('fullpath'));
    if use_old_fslDVARS
        DVARSsh = fullfile(script_dir, 'utilities', 'dvars_nichols.sh');
    else
        dvars_dir = fullfile(script_dir, 'utilities', 'DVARS');
        addpath(dvars_dir);
        addpath(fullfile(dvars_dir, 'Nifti_Util'));
    end
    
    for i = 1:length(subjlist)
        for j = 1:length(tseries)
            tic
            inFile = fullfile(in_dir, subjlist{i}, 'MNINonLinear', 'Results', ...
                tseries{j}, [tseries{j} '.nii.gz']);
            outFile = fullfile(out_dir, subjlist{i}, 'Results', ...
                tseries{j}, [subjlist{i} '_' tseries{j} '_DVARS.txt']);
            disp(['Creating Dvars for ' tseries{i}])
            disp(['Output Directory: ' outFile]); 
            if use_old_fslDVARS
                system([DVARSsh ' ' inFile ' ' outFile]);
            else
                [~, stats] = DVARSCalc(inFile, 'RDVARS');
                try
                    writematrix(stats.RDVARS', outFile);
                catch
                    fid = fopen(outFile, "w");
                    fprintf(fid, '%f\n', stats.RDVARS');
                    fclose(fid);
                end
            end
            toc
        end
    end

    if ~use_old_fslDVARS
        rmpath(dvars_dir);
        rmpath(fullfile(dvars_dir, 'Nifti_Util'));
    end
end