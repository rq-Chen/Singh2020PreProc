function[]=MyStartServerDT_HCP_Simple_Retest(Subjects)
%% Chosen Conditions: WM, Motor, Language, Gambling (b/c gambling is hard)\
if ~iscell(Subjects)
    if ~ischar(Subjects)
        Subjects={num2str(Subjects)};
    end
end

nGSR=3;
copyBase='n';
TaskAll={'WM','MOTOR','LANGUAGE','GAMBLING','RELATIONAL','SOCIAL','EMOTION'};
Atlas='y';
MyStartServerDT_HCP_All_Retest(Subjects,nGSR,copyBase,TaskAll,Atlas);
end