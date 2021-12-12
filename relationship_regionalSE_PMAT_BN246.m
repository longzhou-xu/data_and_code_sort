%% calculte the relationship between SE of ROI and IQ score
clc;clear;close all
clc; clear; close all
%---------------load the MS and SE-------------------------------------
%load the MS and SE of each ROI of each subject
load('D:\PROJECTS\Criticality_cognitive\b_MS_SE_STATIC\MS_SE_BN246_ROI.mat')
Subject_MS_SE = subject;
clear Subject;
%---------------load the HCP cognitive scores---------------------------
load('HCP_iq0804.mat', 'Subject');
Subject_cognitive = Subject;
clear Subject;
load('HCP_iq0804.mat', 'PMAT24_A_CR');
%--------------load the HCP age and edu scores----------------------------
load('subject_age_edu.mat', 'Age_in_Yrs');
age = Age_in_Yrs(2:end);clear Age_in_Yrs;
load('subject_age_edu.mat', 'SSAGA_Educ');
edu = SSAGA_Educ(2:end);clear SSAGA_Educ;
load('subject_age_edu.mat', 'Subject');
Subject_age_edu = Subject(2:end);
clear Subject;
if sum(Subject_age_edu == Subject_cognitive) ~= 1206
    error('wrong');
end
%--------------check and select the subject----------------------------
%-----------demand:(1)same ID; (2)score exits------------------------------
XX =0;
for II = 1 : 295
    for JJ = 1 : 1206
        if Subject_MS_SE(II) == Subject_cognitive(JJ)
            if PMAT24_A_CR(JJ) > 0 && age(JJ) > 0 && edu(JJ) > 0 
                XX = XX+1;
               %------------PMAT------------------------
               PMAT(XX, 1) = II;
               PMAT(XX, 2) = PMAT24_A_CR(JJ);
               PMAT(XX, 3) = age(JJ);
               PMAT(XX, 4) = edu(JJ);
               PMAT(XX, 5) = Subject_MS_SE(II);
               SE_roi(:,XX) = synE(:, II);
            end
        end
    end
end
for ROI=1:246
    plot(SE_roi(ROI,:)',PMAT(:,2),...
        'Marker','o',...
        'MarkerSize',6,...
        'LineStyle','none',...
        'LineWidth',1.5);
    pause(0.2)
    [R_PMAT_SE(ROI,1),P_PMAT_SE(ROI,1)]=corr(SE_roi(ROI,:)',PMAT(:,2));
end
% select the significant ROI, FDR corrected
FDR=mafdr(P_PMAT_SE,'BHFDR', true);
result_PMAT(:,1)=find(FDR<0.05);
result_PMAT(:,2)=R_PMAT_SE(FDR<0.05);
result_PMAT(:,3)=P_PMAT_SE(FDR<0.05);
roi_positive_correct=result_PMAT((result_PMAT(:,2)>0),1);
roi_negative_correct=result_PMAT((result_PMAT(:,2)<0),1);
save('ROI_BN246_PMAT_SE_correlation.mat','R_PMAT_SE','P_PMAT_SE');

% load the BN246.nii to produce the map for BrainNet
addpath('NIfTI_20140122-master')
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