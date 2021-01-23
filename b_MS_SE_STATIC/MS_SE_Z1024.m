% Calculate the MS and SE of Z1024
clc;clear;close all
rest_file=dir('D:\HCP_RS_1_LR\ROI_signal'); % read the folders of subjects
rest_file=rest_file(3:end); % as the first and second are not a effective folders;
for I=1:295
    % read the data
    rest_file(I,1).name
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(I),'.mat']);
    signal=ROI_1024_RS(:,1:1024);
    
    % Z-score normalization
    Cmean(1,:) = mean(signal(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(signal(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
    %calculate the MS and SE
    [syn(I,1),synE(I,1),KOP(:,I)] = syn_synEntropy(signal_z_score(1:1200,:),1200,1024,30);
    subject(I,1)=str2double(rest_file(I,1).name);
end
save('MS_SE_Z1024_RS.mat','subject','syn','synE','KOP');