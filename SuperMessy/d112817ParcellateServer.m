function[]=d112817ParcellateServer(Subs,GSR,type,Parc,varargin)
bDir='/scratch2/Singh/HCP/GSR';
bEnd='/Out/outputs/';
if isempty(Parc)
    Parc='G';
    disp('Assuming Gordon')
end
if strcmpi(Parc(1),'G')
    ParcName='';
else
    ParcName=strcat('_',Parc(1),'_');
end

if isempty(varargin)
Ending='_scrub20_frames400_AveRun1_fix1_pearsonFC_runs_dtseries.mat';
NameEnd='';
else
    NameEnd=strcat('_',varargin{1});
    Ending=strcat('_scrub30_frames150_AveRun1_fix0_pearsonFC_runs_dtseries',varargin{1},'.mat');
end
%Subs=[100307 132017];
if GSR==3
    Mid='_HCP_norm_motion34_GSR';
elseif GSR==2
    Mid='_HCP_norm_motion31_GSR';
elseif GSR==0
    Mid='_HCP_norm_motion25_GSR';
end
for i=1:numel(Subs)
    disp(i)
    load(strcat(bDir,num2str(GSR),bEnd,'sub',num2str(Subs(i)),Mid,num2str(GSR),Ending));
    cc=cell(1,4);
    QC=myQC;
    for j=1:numel(dt)
        Qmod=QC;
        dt{j}=zscore(dt{j}')'; %#ok<AGROW>
        Qmod.run=Qmod.run(j);
        if ~isempty(type)&&strcmpi(type(1),'M')
%        Qmod.run.tmask=qTmask((1+(j-1)*length(qTmask)/4):(j*length(qTmask)/4));    
if strcmpi(Parc(1),'G')
cc{j}=GordonParcellate(dt{j},Qmod,'Ma','Ma','Av','Av','n','y','y');
elseif strcmpi(Parc(1),'M')
cc{j}=MMPParcellate(dt{j},Qmod,'Ma','Ma','Av','Av','n','y','y'); 
else
    error('Parcellations should be G/M')
end
        else
if strcmpi(Parc(1),'G') 
cc{j}=GordonParcellate(dt{j},Qmod,'Av','Av','Av','Av','n','y','y');
elseif strcmpi(Parc(1),'M')
cc{j}=MMPParcellate(dt{j},Qmod,'Av','Av','Av','Av','n','y','y'); 
else
    error('Parcellations should be G/M')
end
        end
    end
    X.Dat=cc;
    X.QC=QC;
    X.Proc={'Av','Av','Av','Av','n','y','y'};
    X.GSR=GSR;
if ~isempty(type)
    save(strcat(bDir,num2str(GSR),'/sub',num2str(Subs(i)),ParcName,type,NameEnd,'.mat'),'X');
else
    save(strcat(bDir,num2str(GSR),'/sub',num2str(Subs(i)),ParcName,NameEnd,'.mat'),'X');
end

end
end
        