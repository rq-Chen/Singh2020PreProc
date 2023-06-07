function[Mask]=GetBrainMask(Atlas)


BaseDir='C:\Users\Matthew\Desktop\DMCC\ATLASES\';
if strcmpi(Atlas(1),'m')
dirMap='HCP-MMP1_on_MNI152_ICBM2009a_nlin_2p4.nii';
    disp(strcat('Getting Mask:',dirMap))
    Map=load_nii(strcat(BaseDir,dirMap));
    Mask=Map.img;
%Ldir='HCP-MMP1_L_on_MNI152_ICBM2009a_nlin_222.nii';
%Rdir='HCP-MMP1_R_on_MNI152_ICBM2009a_nlin_222.nii';

%disp(strcat('Getting Mask:',BaseDir,Ldir))    
%LMap=load_nii(strcat(BaseDir,Ldir));
%disp(strcat('Getting Mask:',BaseDir,Rdir))    
%RMap=load_nii(strcat(BaseDir,Rdir));
%Mask=LMap.img+RMap.img;

disp('Obtained Mask')
elseif strcmpi(Atlas(1),'g')
    BaseDir='C:\Users\Matthew\Desktop\HCP\ATLASES\';
    dirMap='gordon_222_resampled_wsubcort_LPI.nii.gz';
    disp(strcat('Getting Mask:',dirMap))
    Map=load_nii(strcat(BaseDir,'gordon\',dirMap));
    Mask=Map.img;
    Mask(Mask>333)=Mask(Mask>333)-16;
%    Mask=Mask(end:-1:1,:,:);
    disp('Obtained Mask')
else
    error('Atlas should either be MMP or Gordon')
end
end