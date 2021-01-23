%% calculte the relationship between SE of ROI and IQ score
clc;clear;close all
addpath('NIfTI_20140122-master')
%load the MS and SE of each ROI of each subject
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_BN246_ROI.mat'])
subject=subject(1:295,1);
MS_roi=syn(:,1:295);
SE_roi=synE(:,1:295);
load('HCP_iq0804.mat', 'PMAT24_A_CR')
load('HCP_iq0804.mat', 'Subject')
for sub2=1:length(subject)
    for sub1 = 1 : length(Subject)
        if subject(sub2) == Subject(sub1)
            PMAT(sub2,1) = PMAT24_A_CR(sub1);
        end
    end
end
clear Subject PMAT24_A_CR
for ROI=1:246
    plot(SE_roi(ROI,PMAT>-10000)',PMAT(PMAT>-10000),...
        'Marker','o',...
        'MarkerSize',6,...
        'LineStyle','none',...
        'LineWidth',1.5);
    pause(1)
    [R_PMAT_SE(ROI,1),P_PMAT_SE(ROI,1)]=corr(SE_roi(ROI,PMAT>-10000)',PMAT(PMAT>-10000));
end
save('ROI_BN246_PMAT_SE_correlation.mat','R_PMAT_SE','P_PMAT_SE');

% select the significant ROI, FDR corrected
FDR=mafdr(P_PMAT_SE,'BHFDR', true);
result_PMAT(:,1)=find(FDR<0.05);
result_PMAT(:,2)=R_PMAT_SE(FDR<0.05);
result_PMAT(:,3)=P_PMAT_SE(FDR<0.05);
roi_positive_correct=result_PMAT((result_PMAT(:,2)>0),1);
roi_negative_correct=result_PMAT((result_PMAT(:,2)<0),1);

% load the BN246.nii to produce the map for BrainNet
map1=zeros(91,109,91);
map2=zeros(91,109,91);
nii=load_nii(['BN_Atlas_246_2mm.nii']);
mask=nii.img;
for x=1:91
    for y=1:109
        for z=1:91
            if mask(x,y,z) > 0
               map1(x,y,z) = R_PMAT_SE(mask(x,y,z),1);
               for m=1:length(result_PMAT)
                   if mask(x,y,z) == result_PMAT(m,1)
                       map2(x,y,z) = R_PMAT_SE(mask(x,y,z),1);
                   end
               end
            end
        end
    end
end
nii.img=map1;
save('relationship_regional_SE_PMAT_BN246.mat', 'map1');
nii.hdr.dime.bitpix=32;
nii.hdr.dime.datatype=16;
save_nii(nii, 'ROI_R_PMAT_SE_BN246.nii')
nii.img=map2;
save_nii(nii, 'ROI_R_PMAT_SE_BN246_significant.nii')