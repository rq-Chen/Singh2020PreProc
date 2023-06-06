function[]=NewDVARScalc(InDir,OutDir)
%% Use to replace old method of calculating DVARS (shell) b/c won't run on server
A=DVARScalc(InDir,'RDVARS','Norm',1);
A=A-min(A)+1;
dlmwrite(OutDir,A);
end
