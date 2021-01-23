% calculate the MS and SE in each ROI using the vox-wise data
clc; clear; close all;
% read the data folder of each subject
file_dir=dir('E:\fMRI_signal\HCP_REST_LR\BN_246_VOX_signal');
file_dir=file_dir(3:end); % as the first and second clumun are not effect folders
% calculate the MS and SE of each ROI of each subject
for Sub=1:295
    for ROI=1:246
        load(['E:\fMRI_signal\HCP_REST_LR\BN_246_VOX_signal\',...
            file_dir(Sub).name,'\ROI',num2str(ROI),...
            '\RS_BN_246.mat']);
        signal=RS_BN_246';
    
        % Z-score normalization
        Cmean(1,:) = mean(signal(1:1200,:),1);
        Cmeanmatrix = repmat(Cmean,1200,1);
        Cstd = std(signal(1:1200,:),1);
        Cstdmatrix = repmat(Cstd,1200,1);
        signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
        %calculate the MS and SE
        cell_number=length(signal(1,:));
        [syn(ROI,Sub),synE(ROI,Sub)] = syn_synEntropy(signal_z_score(1:1200,:),1200,cell_number,30);
        subject(Sub,1)=str2double(file_dir(Sub,1).name);
        clear signal Cmean Cmeanmatrix Cstd Cstdmatrix signal_z_score
    end
end
save('MS_SE_BN246_ROI.mat','syn','synE','subject');
%%
clc;clear;close all;
% read the data folder of each subject
file_dir=dir('E:\fMRI_signal\HCP_REST_LR\HY_96_VOX_signal');
file_dir=file_dir(3:end); % as the first and second clumun are not effect folders
% calculate the MS and SE of each ROI of each subject
for Sub=1:295
    for ROI=1:96
        Sub
        file_dir(Sub).name
        ROI
        load(['E:\fMRI_signal\HCP_REST_LR\HY_96_VOX_signal\',...
            file_dir(Sub).name,'\ROI',num2str(ROI),...
            '\RS_HY_96.mat']);
        signal=RS_HY_96';
    
        % Z-score normalization
        Cmean(1,:) = mean(signal(1:1200,:),1);
        Cmeanmatrix = repmat(Cmean,1200,1);
        Cstd = std(signal(1:1200,:),1);
        Cstdmatrix = repmat(Cstd,1200,1);
        signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
        %calculate the MS and SE
        cell_number=length(signal(1,:));
        [syn(ROI,Sub),synE(ROI,Sub)] = syn_synEntropy(signal_z_score(1:1200,:),1200,cell_number,30);
        subject(Sub,1)=str2double(file_dir(Sub,1).name);
        clear signal Cmean Cmeanmatrix Cstd Cstdmatrix signal_z_score
    end
end
save('MS_SE_HY96_ROI.mat','syn','synE','subject');
%%
clc;clear;close all
% read the data folder of each subject
file_dir=dir('E:\fMRI_signal\HCP_REST_LR\Z_1024_VOX_signal');
file_dir=file_dir(3:end); % as the first and second clumun are not effect folders
% calculate the MS and SE of each ROI of each subject
for Sub=1:295
    for ROI=1:1024
        Sub
        file_dir(Sub).name
        ROI
        load(['E:\fMRI_signal\HCP_REST_LR\Z_1024_VOX_signal\',...
            file_dir(Sub).name,'\ROI',num2str(ROI),...
            '\RS_Z_1024.mat']);
        signal=RS_Z_1024';
    
        % Z-score normalization
        Cmean(1,:) = mean(signal(1:1200,:),1);
        Cmeanmatrix = repmat(Cmean,1200,1);
        Cstd = std(signal(1:1200,:),1);
        Cstdmatrix = repmat(Cstd,1200,1);
        signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
        %calculate the MS and SE
        cell_number=length(signal(1,:));
        [syn(ROI,Sub),synE(ROI,Sub)] = syn_synEntropy(signal_z_score(1:1200,:),1200,cell_number,30);
        subject(Sub,1)=str2double(file_dir(Sub,1).name);
        clear signal Cmean Cmeanmatrix Cstd Cstdmatrix signal_z_score
    end
end
save('MS_SE_Z1024_ROI.mat','syn','synE','subject');