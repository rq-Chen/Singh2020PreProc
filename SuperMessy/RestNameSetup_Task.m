function tseries = RestNameSetup_Task(Task)

    if strcmpi(Task(1),'W')||strcmpi(Task(1),'N')
       tseries = {'tfMRI_WM_LR', 'tfMRI_WM_RL'};
    elseif strcmpi(Task(1),'M')
        tseries = {'tfMRI_MOTOR_LR', 'tfMRI_MOTOR_RL'};
    else
        tseries = {strcat('tfMRI_',Task,'_LR'),strcat('tfMRI_',Task,'_RL')};                  
    end
%    if strcmp( DataType, 'HCP' )%if it's HCP Data Setup these folders
%        tseries = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL',...
%                   'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}
%    else%if the data is DMCC 
%          tseries = {'rfMRI_RestBas1_AP', 'rfMRI_RestBas2_PA', 'rfMRI_RestPro1_AP', ...
%            'rfMRI_RestPro2_PA' , 'rfMRI_RestRea1_AP','rfMRI_RestRea2_PA'}
%    end
    
end