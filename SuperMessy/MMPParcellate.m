function[Out]=MMPParcellate(X,QC,Corttype,Subtype,Stemtype,Ceretype,doNorm,doInterp,isServer,varargin)

%% type arguments for parcellating: 'Av','Ma','PC'
%% Ave computes average
%% Mahal does Mahalanobis Distance
%% PC1 does first principal component
% Parc.NP = 360;%Number of Parcels
% Parc.name = 'MMP360';%Name of the parcel will be used to store the ptseries
% Parc.ordered = 0;%is the parcel already ordered Correctly ie. Community then Hemisphere
% Parc.key = '/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/HCP-MMP/Glasser_et_al_2016_HCP_MMP1.0_RVVG/MMP360ParcelsKey.csv'

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


Combo=[gCort;MatchSub(G,SubCifti,ParcInfo)];

if strcmpi(doNorm(1),'y')
    X=zscore(X')';
end



%whos Combo
%whos X
Out=zeros(max(Combo),size(X,2));
%unique(Combo(Combo~=0))
for i=unique(Combo(Combo~=0)')
    disp(i)
    if i>ParcInfo.NP
        if or(i==(ParcInfo.NP+16),i==(ParcInfo.NP+17))
            type=Ceretype;
        elseif i==(ParcInfo.NP+15)
            type=Stemtype;
        else
        type=Subtype;
        end
    else
        type=Corttype;
    end
    if strcmpi(type(1:2),'Av')
        Out(i,:)=mean(X(Combo==i,:),1);
    elseif strcmpi(type(1:2),'Ma')
        Cov=cov(X(Combo==i,:)');
        Out(i,:)=QuadMat(X(Combo==i,:),pinv(Cov));
    elseif strcmpi(type(1:2),'PC')
        Cov=cov(X(Combo==i,:)');
        [u,~]=eig(Cov);
        Out(i,:)=u(:,end)'*X(Combo==i,:);
    else
        error('Type should be Av Ma or PC')
    end
end
%% Removed this line on 5.17.18--need to redo all files to match!
%Out(1:ParcInfo.NP,:)=reorganizeSRI(Out(1:ParcInfo.NP,:),ParcInfo);
        
if strcmpi(doInterp(1),'y')
    tt=[QC.run.tmask];
    Out(:,tt==0)=Out;
    Out(:,tt==1)=interp1(find(tt==0),Out(:,tt==0)',find(tt==1))';
end
%% Remove first frame
Out(:,1)=[];

%% Used for matching subcorticals in Mask to the volume
%% So Doesn't matter if using the Gordon Mask
function[fOut]=MatchSub(VolMask,CiftiMask,ppInf)
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


end
