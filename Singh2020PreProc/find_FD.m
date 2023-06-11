% 2023/06/09 Modified by Ruiqi Chen:
%
% Now the function only searches for FD in the input and NOT in
% the output. If FD is not found in the input, it will generates
% it and overwrite any existing one in the output.
%
% --------------------------------------------------------------
%
%%This funcion Checks for the FD's in both the input and the output.
%if it is in the input copy it to the output
%if the FD cannot be found generate them for the output
function find_FD(in, out, Subjects, tseries)
    if ~iscell(Subjects)
        Subjects = {Subjects};
    end
    for i = 1:numel(Subjects)
        for j = 1:numel(tseries)
        
            Subject = Subjects{i};
            input = fullfile(in, Subject,'MNINonLinear','Results', tseries{j});%the file where the FDs should be
            output = fullfile(out, Subject, 'Results', tseries{j} );%the file where we want to copt the FDs to
            FDName = [Subject '_' tseries{j} '_FD.txt']%name for the FDs
            if exist(fullfile(input, FDName), 'file')%is the file in the input
                disp(['Found ' FDName ' in the input: ' input ])
                copyfile(fullfile(input, FDName), fullfile(output, FDName));%copy it to the output

            % elseif exist(fullfile(output,FDName), 'file')%is the file in the output
            %     disp(['Found ' FDName ' in the output: ' output])
            
            else 
                disp(['Cannot find ' FDName ' in the input. ']) %generate the FD's
                % disp([output ' or input'])
                % disp(input)
                disp('Will Begin generating')
                disp(['Placing in: ' output])
                c_calc_FD(in, output, {Subject}, tseries(j))
            
            end
        end
    end
end
    
                
        