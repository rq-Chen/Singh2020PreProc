% bandstop_filter
% for [0.06 - 0.14] Hz, removing respiratory frequencies (Siegel)
% used for DMCC Data
% 
function filter_FD_DVARS(Filter, Subject, in_dir, out_dir, tseries)

% FDs

    sub_id = num2str(Subject)
    for i=1:length(tseries)
        %Read in the FD and Filter
        out = fullfile(out_dir, sub_id, 'Results', tseries{i} );%the file where we want to copt the FDs to
        if ~isfile(fullfile(out, [sub_id '_' tseries{i} '_FD_FILT.txt']))
            FD = [sub_id '_' tseries{i} '_FD.txt']
            if ~isfile(fullfile(out, FD))
                warning([tseries{i} ' FD not found!'])
                continue
            end
            FD_ = dlmread(fullfile(out, FD));%read in the FD
            FD_FILT = filtfilt(Filter,FD_);%filter the FD
            dlmwrite(fullfile(out, [sub_id '_' tseries{i} '_FD_FILT.txt']), FD_FILT)%write the filtered output
        end
        
        % Copy the DVARs from the input Directory to the output Directory
        % Read in the DVARS and Filter repeat the same thing as above with the
        % DVars
        if ~isfile(fullfile(out, [sub_id '_' tseries{i} '_DVARS_FILT.txt']))
            DVAR = [sub_id '_' tseries{i} '_DVARS.txt']
            if ~isfile(fullfile(out, DVAR))
                warning([tseries{i} ' DVARS not found!'])
                continue
            end
            DVARS = dlmread(fullfile(out,DVAR));
            DVARS_FILT = filtfilt(Filter,DVARS);    
            dlmwrite(fullfile(out, [sub_id '_' tseries{i} '_DVARS_FILT.txt']), DVARS_FILT)
        end
    end
end