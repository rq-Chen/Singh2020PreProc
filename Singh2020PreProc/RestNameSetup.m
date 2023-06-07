function tseries = RestNameSetup(DataType)
    if strcmp( DataType, 'HCP' )%if it's HCP Data Setup these folders
        tseries = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL',...
                   'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}
    else%if the data is DMCC 
        tseries = {'rfMRI_RestBas1_AP', 'rfMRI_RestBas2_PA', 'rfMRI_RestPro1_AP', ...
            'rfMRI_RestPro2_PA', 'rfMRI_RestRea1_AP','rfMRI_RestRea2_PA' }
    end
end