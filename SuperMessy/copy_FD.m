%Copy the FD's to the output directory
% 
function copy_FD(Subjlist, in_dir, out_dir, tseries)

    for i = 1:length(Subjlist)

        sub_id = num2str(Subject)
        for i=1:length(tseries)
            %Copy the FDs from the input Directory to the Output Directory
            %Read in the FD and Filter
            in = fullfile(in_dir, sub_id,'MNINonLinear','Results', tseries{i});%the file where the FDs should be
            out = fullfile(out_dir, sub_id, 'Results', tseries{i} );%the file where we want to copt the FDs to
            FD = [sub_id '_' tseries{i} '_FD.txt']%name for the FDs
    
            try
                copyfile(fullfile(in, FD),fullfile(out, FD)); % copy the file
            catch
                error(['Could Not Filter FDs, Check to make sure that they are present in the following folder: ' fullfile(in, FD)])
            end
        end
    end
end

