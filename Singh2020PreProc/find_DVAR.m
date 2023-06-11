% 2023/06/09 Modified by Ruiqi Chen:
%
% Now the function only searches for DVARS in the input and NOT in
% the output. If DVARS is not found in the input, it will generates
% it and overwrite any existing one in the output.
%
% --------------------------------------------------------------
%
%This funcion Checks for the DVAR's in both the input and the output.
%if it is in the input copy it to the output
%if the DVARs cannot be found generate them for the output
function find_DVAR(in, out, Subject, tseries)
    for i = 1:length(Subject)
        for j = 1:length(tseries)
            
            input = fullfile(in, Subject,'MNINonLinear','Results', tseries{j});%the file where the DVARs should be
            output = fullfile(out, Subject, 'Results', tseries{j} );%the file where we want to copt the DVARs to
            DVARName = [Subject '_' tseries{j} '_DVARS.txt']%name for the DVARs
                    
            if exist(fullfile(input, DVARName), 'file')%is the file in the input
                disp(['Found ' DVARName ' in the input: ' input ])
                copyfile(fullfile(input, DVARName), fullfile(output, DVARName));%copy it to the output

            else
                % if exist(fullfile(output,DVARName), 'file');             %is the file in the output
                %     file=dir(fullfile(output,DVARName));
                %     dvSize=file.bytes;
                % else
                %     dvSize=0;
                % end
                % if dvSize>100

                %     disp(['Found ' DVARName ' in the output: ' output])
                
                % else
                    
                    disp(['Cannot find ' DVARName ' in the input. Will begin Generating...']) %generate the DVAR's
                    % disp([output ' or input'])
                    % disp(input)
                    % disp('Will Begin generating')
                    disp(['Placing in: ' output])
                    
                    calc_DVAR_mod(in, out, {Subject}, tseries(j))
                % end
            end
        end
    end