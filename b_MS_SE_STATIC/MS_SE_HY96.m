clc;clear;close all
% Calculate the MS and SE of HY96
for S = 1:295
    % read the data
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S),'.mat']);
    signal=rest_HY96_ROI(:,1:96);
    % Z-score normalization
    Cmean(1,:) = mean(signal(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(signal(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    %calculate the MS and SE
    [syn(S,1),synE(S,1),KOP(:,S)] = syn_synEntropy(signal_z_score(1:1200,:),1200,96,30);
end
save('MS_SE_HY96_RS.mat','syn','synE','KOP');