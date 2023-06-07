function[Out]=GordonParcellate(X,QC,Corttype,Subtype,Stemtype,Ceretype,doNorm,doInterp,isServer,varargin)

%% type arguments for parcellating: 'Av','Ma','PC'
%% Ave computes average
%% Mahal does Mahalanobis Distance
%% PC1 does first principal component
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


Combo=[gCort;MatchSub(G,SubCifti)];

if strcmpi(doNorm(1),'y')
    X=zscore(X')';
end



%whos Combo
%whos X
Out=zeros(max(Combo),size(X,2));
%unique(Combo(Combo~=0))
for i=unique(Combo(Combo~=0)')
    disp(i)
    if i>333
        if or(i==349,i==350)
            type=Ceretype;
        elseif i==348
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
Out(1:333,:)=reorganizeSRI(Out(1:333,:),ParcInfo);
        
if strcmpi(doInterp(1),'y')
    tt=[QC.run.tmask];
    Out(:,tt==0)=Out;
    Out(:,tt==1)=interp1(find(tt==0),Out(:,tt==0)',find(tt==1))';
end
%% Remove first frame
Out(:,1)=[];


function[fOut]=MatchSub(VolMask,CiftiMask)
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


end
