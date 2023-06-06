function[Out]=GetAtlasMasks(Atlas,isServer)
%% Get Masks for atlases

if strcmpi(Atlas(1),'G')
    
ParcInfo.NP=333;
if strcmpi(isServer(1),'n')
G=GetBrainMask('g');
ParcInfo.key='C:\Users\Matthew\Desktop\HCP\Fixed Atlas\Gordon333SRI_ParcelsKey.csv';
SubCifti=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\HCP_subcortical\Subcortnew.txt');
gCort=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\MyMask\newGord.txt');
else
    if strcmpi(isServer(1),'r')
        mDir='/home/SGE/f.singh/ATLASES/MyMask/';
    else
    mDir='/scratch2/Singh/ATLASES/MyMask/';
    end
    G=load(strcat(mDir,'GordonWsub.mat'));
    G=G.G;
    ParcInfo.key=strcat(mDir,'Fixed Atlas/Gordon333SRI_ParcelsKey.csv');
    gCort=importdata(strcat(mDir,'MyMask/newGord.txt'));
    SubCifti=importdata(strcat(mDir,'HCP_subcortical/Subcortnew.txt'));
end


Combo=[gCort;MatchSubG(G,SubCifti)];
Out={Combo,ParcInfo};
elseif strcmpi(Atlas(1),'M')
    
ParcInfo.NP=360;
if strcmpi(isServer(1),'n')
G=GetBrainMask('g');
%ParcInfo.key='C:\Users\Matthew\Desktop\HCP\Fixed Atlas\Gordon333SRI_ParcelsKey.csv';

SubCifti=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\HCP_subcortical\Subcortnew.txt');
gCort=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\MyMask\newMMP.txt');
else
    if strcmpi(isServer(1),'r')
        mDir='/home/SGE/f.singh/ATLASES/MyMask/';
    else
    mDir='/scratch2/Singh/ATLASES/MyMask/';
    end
    G=load(strcat(mDir,'GordonWsub.mat')); %% Doesn't matter that using gordon here--only for the subcorticals anyway
    G=G.G;
    G(G>0)=G(G>0)+(ParcInfo.NP-333);
    ParcInfo.key=strcat(mDir,'Fixed Atlas/MMP360_ParcelsKey.csv');
    gCort=importdata(strcat(mDir,'MyMask/MMPlab.txt'));
    SubCifti=importdata(strcat(mDir,'HCP_subcortical/Subcortnew.txt'));
end


Combo=[gCort;MatchSubM(G,SubCifti,ParcInfo)];
Out={Combo,ParcInfo};
elseif strcmpi(Atlas(1),'Y')
YeoSize={'01','02','04','06','08','10'};
ParcName=cell(1,numel(YeoSize));
Out=cell(1,numel(YeoSize));
    for iY=1:numel(YeoSize)
        ParcName{iY}=strcat('Y',YeoSize{iY});
        nYeo=strcat('Yeo_',YeoSize{iY});
        ParcInfo.NP=100*str2double(nYeo([end-1 end]));
if strcmpi(isServer(1),'n')
%G=GetBrainMask('g');
G=load('C:\Users\Matthew\Documents\MATLAB\GordonWsub.mat');
    G=G.G;
    G(G>0)=G(G>0)+(ParcInfo.NP-333);
%ParcInfo.key='C:\Users\Matthew\Desktop\HCP\Fixed Atlas\Gordon333SRI_ParcelsKey.csv';

SubCifti=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\HCP_subcortical\Subcortnew.txt');
gCort=importdata(strcat('C:\Users\Matthew\Desktop\HCP\ATLASES\Yeo\',nYeo,'.txt'));
load('C:\Users\Matthew\Desktop\HCP\ATLASES\Yeo\cortInd.mat','cortInd');
else
    if strcmpi(isServer(1),'r')
        mDir='/home/SGE/f.singh/ATLASES/MyMask/';
    else
    mDir='/scratch2/Singh/ATLASES/Yeo/';
    end
    G=load(strcat(mDir,'GordonWsub.mat')); %% Doesn't matter that using gordon here--only for the subcorticals anyway
    load(strcat(mDir,'cortInd.mat'),'cortInd');
    G=G.G;
    G(G>0)=G(G>0)+(ParcInfo.NP-333);
   % ParcInfo.key=strcat(mDir,'Fixed Atlas/MMP360_ParcelsKey.csv');
    gCort=importdata(strcat(mDir,nYeo,'.txt'));
    SubCifti=importdata(strcat(mDir,'HCP_subcortical/Subcortnew.txt'));
end
gCort=gCort(cortInd);
Out{iY}=[gCort;MatchSubY(G,SubCifti,ParcInfo)];
    end
end

function[fOut]=MatchSubG(VolMask,CiftiMask)
uu=unique(VolMask);
uu=uu(uu>333);
VoxCount=zeros(1,numel(uu));
for fi=1:numel(uu)
    VoxCount(fi)=sum(sum(sum(VolMask==uu(fi))));
end
fOut=zeros(size(CiftiMask));
cc=unique(CiftiMask);
for fi=1:numel(cc)
    [~,kk]=min(abs(VoxCount-sum(CiftiMask==cc(fi))));
    fOut(CiftiMask==cc(fi))=uu(kk);
end
end

function[fOut]=MatchSubM(VolMask,CiftiMask,ppInf)
uu=unique(VolMask);
uu=uu(uu>ppInf.NP);
VoxCount=zeros(1,numel(uu));
for fi=1:numel(uu)
    VoxCount(fi)=sum(sum(sum(VolMask==uu(fi))));
end
fOut=zeros(size(CiftiMask));
cc=unique(CiftiMask);
for fi=1:numel(cc)
    [~,kk]=min(abs(VoxCount-sum(CiftiMask==cc(fi))));
    fOut(CiftiMask==cc(fi))=uu(kk);
end
end

function[fOut]=MatchSubY(VolMask,CiftiMask,ppInf)
uu=unique(VolMask);
uu=uu(uu>ppInf.NP);
VoxCount=zeros(1,numel(uu));
for fi=1:numel(uu)
    VoxCount(fi)=sum(sum(sum(VolMask==uu(fi))));
end
fOut=zeros(size(CiftiMask));
cc0=unique(CiftiMask);
for fi=1:numel(cc0)
    [~,kk]=min(abs(VoxCount-sum(CiftiMask==cc0(fi))));
    fOut(CiftiMask==cc0(fi))=uu(kk);
end
end

end