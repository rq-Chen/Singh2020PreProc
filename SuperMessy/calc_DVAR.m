function calc_DVAR(in_dir, out_dir, subjlist, tseries)
    for i = 1:length(subjlist)
        for j = 1:length(tseries)
            tic
            output =  fullfile(out_dir, subjlist{i},'Results', tseries{j}, [subjlist{i} '_' tseries{j} '_DVARS.txt']);
            disp(['Creating Dvars for ' tseries{i}])
            disp(['Output Directory: ' output]); 
            system(['/scratch1/MitchJeffers/restingState/dvars_nichols.sh ' fullfile(in_dir, subjlist{i}, ...
                '/MNINonLinear/Results',tseries{j}, [tseries{j} '.nii.gz ']) output]);              
            toc
        end
    end
end