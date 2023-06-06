%This funcion Checks for the FD's in both the input and the output.
%if it is in the input copy it to the output
%if the FD cannot be found generate them for the output
function find_FDretest(in, out, Subject, tseries)
        for j = 1:length(tseries)
            
            input = fullfile(in, strcat(Subject,'_3T/RESOURCES/',tseries{j},'_PostFix'),'MNINonLinear','Results', tseries{j});%the file where the FDs should be
            output = fullfile(out, Subject, 'Results', tseries{j} );%the file where we want to copt the FDs to
            FDName = [Subject '_' tseries{j} '_FD.txt']%name for the FDs
                    
            if exist(fullfile(input, FDName), 'file')%is the file in the input
                disp(['Found ' FDName ' in the input: ' input ])
                copyfile(fullfile(input, FDName), fullfile(output, FDName));%copy it to the output

            elseif exist(fullfile(output,FDName), 'file')%is the file in the output
                disp(['Found ' FDName ' in the output: ' output])
            
            else 
                disp(['Cannot find ' FDName ' in either:']) %generate the FD's
                disp([output ' or input'])
                disp(input)
                disp('Will Begin generating')
                disp(['Placing in: ' output])
                c_calc_FDretest(in, output, {Subject}, {tseries{j}})
            
            end
        end
    end
    
                
        