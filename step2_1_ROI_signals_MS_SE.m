%% step2: calculate the mean synchronization (MS) and synchronization entropy(SE)
% setting environment and workspace
clc;clear;close all
addpath('function')

%% HY96-ROI signals
% setting environment and workspace
clc;clear;close all

% calculate the MS, SE and Kuramoto order paramater
for SUB = 1:295
    % the ROI signal path
    HY96_ROIs_path = ['step_1\ROI_signals\HY_96\sub',num2str(SUB),'.mat'];
    % load the ROI-signals
    load(HY96_ROIs_path);
    signals = rest_HY96_ROI(:,1:96);% the ROI signals
    clear rest_HY96_ROI
    
    % z-score normalized signals
    signals_zs = zscore(signals);
    %calculate the MS, SE and kuramoto parameter
    time_len = 1200;
    node_num = 96;
    entropy_bin = 30;
    [MS(SUB,1),SE(SUB,1),KOP(:,SUB)] = syn_synEntropy(signals_zs, time_len, node_num, entropy_bin);
end
save('MS_SE_HY96_RS.mat','MS','SE','KOP');
%% BN246-ROI signals
% setting
clc;clear;close all
% calculate the MA, SE and kuramoto order parameter
for SUB=1:295
    % the BN246 ROI signals path
    clear BN246_ROIs_path
    BN246_ROIs_path = ['step_1\ROI_signals\BA_246\sub',num2str(SUB),'.mat'];
    % load the BN246 signals
    load(BN246_ROIs_path);
    signals = ROI_246_RS;
    clear ROI_246_RS
    
    % Z-score normalization
    signals_zs = zscore(signals);
    %calculate the MS and SE
    time_len = 1200;
    node_num = 246;
    entropy_bin = 30;
    [MS(SUB,1), SE(SUB,1), KOP(:,SUB)] = syn_synEntropy(signals_zs,time_len,node_num,entropy_bin);
end
save('MS_SE_BA246_RS.mat','MS','SE','KOP');
%% BN210-ROI signals
% setting
clc;clear;close all
% calculating
for SUB=1:295
    % set the data path
    BN246_ROIs_path = ['step_1\ROI_signals\BA_246\sub',num2str(SUB),'.mat'];
    % load the data
    load(BN246_ROIs_path);
    signals = ROI_246_RS(:, 1:210);
    clear ROI_246_RS
    
    % Z-score normalization
    signals_zs = zscore(signals);
    %calculate the MS and SE
    time_len = 1200;
    node_num = 210;
    entropy_bin = 30;
    [MS(SUB,1),SE(SUB,1),KOP(:,SUB)] = syn_synEntropy(signals_zs, time_len, node_num, entropy_bin);
end
save('MS_SE_BA210_RS.mat','MS','SE','KOP');

%% Z-1024 ROI signals
% setting
clc;clear;close all
% calculating
for SUB=1:295
    % set the data path
    Z1024_ROIs_path = ['step_1\ROI_signals\Z_1024\sub',num2str(SUB),'.mat']; 
    load(Z1024_ROIs_path);
    signals=ROI_1024_RS(:,1:1024);
    
    % Z-score normalization
    signals_zs = zscore(signals);
    
    %calculate the MS and SE
    time_len =1200;
    node_num = 1024;
    entropy_bin = 30;
    [MS(SUB,1),SE(SUB,1),KOP(:,SUB)] = syn_synEntropy(signals_zs, time_len, node_num, entropy_bin);
end
save('MS_SE_Z1024_RS.mat','MS','SE','KOP');

