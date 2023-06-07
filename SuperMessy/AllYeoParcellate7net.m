function [] = AllYeoParcellate7net(Subs,GSR,isServer,type,doMSM,varargin)
%% Performs all levels of Yeo parcellation
doInterp='y';
YeoSize={'01','02','04','06','08','10'};
Corttype='Av';Subtype='Av';Stemtype='Av';Ceretype='Av';
if strcmpi(isServer(1),'n')
    bDir='C:\Users\Matthew\Desktop\HCP\Siegel\GSR';
else
if strcmpi(doMSM(1),'y')
bDir='/scratch2/Singh/HCP/SuperClean/GSR';
else    
bDir='/scratch2/Singh/HCP/GSR';
end
end
bEnd='/Out/outputs/';

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

for iSub=1:numel(Subs)
    disp(iSub)
    load(strcat(bDir,num2str(GSR),bEnd,'sub',num2str(Subs(iSub)),Mid,num2str(GSR),Ending),'myQC','QC','Subjlist','switches','excluded','dt');
    for iY=1:numel(YeoSize)
        ParcName=strcat('Y',YeoSize{iY});
    cc=cell(1,4);
    QC=myQC;
    for j=1:numel(dt)
        dt{j}=zscore(dt{j}')'; %#ok<AGROW>
        Qmod=QC;
%        dt{j}=zscore(dt{j}')'; %#ok<AGROW> %% Already normalize elsewhere
        Qmod.run=Qmod.run(j);
nYeo=strcat('Yeo_',YeoSize{iY});       
%nYeo=strcat(nYeo(1:3),'_',nYeo([end-1 end]));
ParcInfo.NP=100*str2double(nYeo([end-1 end]));
%ParcInfo.NP=360;
if strcmpi(isServer(1),'n')
%G=GetBrainMask('g');
G=load('C:\Users\Matthew\Documents\MATLAB\GordonWsub.mat');
    G=G.G;
    G(G>0)=G(G>0)+(ParcInfo.NP-333);
%ParcInfo.key='C:\Users\Matthew\Desktop\HCP\Fixed Atlas\Gordon333SRI_ParcelsKey.csv';

SubCifti=importdata('C:\Users\Matthew\Desktop\HCP\ATLASES\HCP_subcortical\Subcortnew.txt');
gCort=importdata(strcat('C:\Users\Matthew\Desktop\DMCC\ATLASES\Yeo_7net\',nYeo,'.txt'));
load('C:\Users\Matthew\Desktop\HCP\ATLASES\Yeo\cortInd.mat','cortInd');
else
    if strcmpi(isServer(1),'r')
        mDir='/home/SGE/f.singh/ATLASES/MyMask/';
    else
    mDir='/scratch2/Singh/ATLASES/Yeo_7net/';
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
Combo=[gCort;MatchSub(G,SubCifti,ParcInfo)];
size(Combo)
%unique(Combo);
Out=zeros(max(Combo),size(dt{j},2));
for i=unique(Combo(Combo~=0)')
    disp(i)
        Out(i,:)=mean(dt{j}(Combo==i,:),1);
end
type=[];
if strcmpi(doInterp(1),'y')
    tt=[Qmod.run.tmask];
    Out(:,tt==0)=Out;
    Out(:,tt==1)=interp1(find(tt==0),Out(:,tt==0)',find(tt==1))';
end
%% Remove first frame
Out(:,1)=[];
cc{j}=Out;
    end
    X.Dat=cc;
    X.QC=QC;
    X.Proc={'Av','Av','Av','Av','n','y','y'};
    X.GSR=GSR;
    if strcmpi(doMSM(1),'y')
    X.Misc='SuperClean';
    end
    X.Nets='Yeo 7 nets';
if ~isempty(type)
    save(strcat(bDir,num2str(GSR),'/sub',num2str(Subs(iSub)),ParcName,type,NameEnd,'_7net.mat'),'X');
else
    save(strcat(bDir,num2str(GSR),'/sub',num2str(Subs(iSub)),ParcName,NameEnd,'_7net.mat'),'X');
end

    end
end

function[fOut]=MatchSub(VolMask,CiftiMask,ppInf)
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

