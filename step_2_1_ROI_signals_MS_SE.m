%% step2: calculate the mean synchronization (MS) and synchronization entropy(SE)

%% HY96-ROI signals
clc;
clear;
close all
addpath Function;

% load the ROI signals 

% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'), 'ROIsignals_HY');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 96;
    signals = ROIsignals_HY(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = xlz_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  xlz_kop2sta(kop, bins_num);
    
    kop_sta_HY(sub).subj_ID = ROIsignals_HY(sub).subj_ID;
    kop_sta_HY(sub).kop = kop;
    kop_sta_HY(sub).mean_kop = MS(sub);
    kop_sta_HY(sub).min_kop = min(kop);
    kop_sta_HY(sub).max_kop = max(kop);
    kop_sta_HY(sub).std_kop = SS(sub);
    kop_sta_HY(sub).cv_kop = CS(sub);
    kop_sta_HY(sub).bins_num = bins_num; 
    kop_sta_HY(sub).entropy_kop = SE(sub);
    kop_sta_HY(sub).sample_failed = sample_failed;
    
end

%% BN246-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'),  'ROIsignals_BN');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 246;
    signals = ROIsignals_BN(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = xlz_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  xlz_kop2sta(kop, bins_num);
    
    kop_sta_BN246(sub).subj_ID = ROIsignals_BN(sub).subj_ID;
    kop_sta_BN246(sub).kop = kop;
    kop_sta_BN246(sub).mean_kop = MS(sub);
    kop_sta_BN246(sub).min_kop = min(kop);
    kop_sta_BN246(sub).max_kop = max(kop);
    kop_sta_BN246(sub).std_kop = SS(sub);
    kop_sta_BN246(sub).cv_kop = CS(sub);
    kop_sta_BN246(sub).bins_num = bins_num; 
    kop_sta_BN246(sub).entropy_kop = SE(sub);
    kop_sta_BN246(sub).sample_failed = sample_failed;
    
end

%% BN210-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'),  'ROIsignals_BN');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 210;
    signals = ROIsignals_BN(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = xlz_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  xlz_kop2sta(kop, bins_num);
    
    kop_sta_BN210(sub).subj_ID = ROIsignals_BN(sub).subj_ID;
    kop_sta_BN210(sub).kop = kop;
    kop_sta_BN210(sub).mean_kop = MS(sub);
    kop_sta_BN210(sub).min_kop = min(kop);
    kop_sta_BN210(sub).max_kop = max(kop);
    kop_sta_BN210(sub).std_kop = SS(sub);
    kop_sta_BN210(sub).cv_kop = CS(sub);
    kop_sta_BN210(sub).bins_num = bins_num; 
    kop_sta_BN210(sub).entropy_kop = SE(sub);
    kop_sta_BN210(sub).sample_failed = sample_failed;
    
end

%% Z1024-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'),  'ROIsignals_Z');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 1024;
    signals = ROIsignals_Z(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = xlz_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  xlz_kop2sta(kop, bins_num);
    
    kop_sta_Z(sub).subj_ID = ROIsignals_Z(sub).subj_ID;
    kop_sta_Z(sub).kop = kop;
    kop_sta_Z(sub).mean_kop = MS(sub);
    kop_sta_Z(sub).min_kop = min(kop);
    kop_sta_Z(sub).max_kop = max(kop);
    kop_sta_Z(sub).std_kop = SS(sub);
    kop_sta_Z(sub).cv_kop = CS(sub);
    kop_sta_Z(sub).bins_num = bins_num; 
    kop_sta_Z(sub).entropy_kop = SE(sub);
    kop_sta_Z(sub).sample_failed = sample_failed;
    
end

save(fullfile('step_2_MS_SE_relationship', 'kuramoto_order_paramter.mat'), 'kop_sta_HY', 'kop_sta_BN246', 'kop_sta_BN210', 'kop_sta_Z', '-v7.3')
