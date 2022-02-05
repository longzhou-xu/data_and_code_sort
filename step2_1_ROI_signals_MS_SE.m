%% step2: calculate the mean synchronization (MS) and synchronization entropy(SE)

%% HY96-ROI signals
clc;
clear;
close all
addpath Function;

% load the ROI signals 

% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_HY');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 96;
    signals = ROIsignals_HY(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = f_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  f_kop2sta(kop, bins_num);
    
    results_HY(sub).subj_ID = ROIsignals_HY(sub).subj_ID;
    results_HY(sub).kop = kop;
    results_HY(sub).mean_kop = MS(sub);
    results_HY(sub).min_kop = min(kop);
    results_HY(sub).max_kop = max(kop);
    results_HY(sub).std_kop = SS(sub);
    results_HY(sub).cv_kop = CS(sub);
    results_HY(sub).bins_num = bins_num; 
    results_HY(sub).entropy_kop = SE(sub);
    results_HY(sub).sample_failed = sample_failed;
    
end

%% BN246-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_BN');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 246;
    signals = ROIsignals_BN(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = f_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  f_kop2sta(kop, bins_num);
    
    results_BN246(sub).subj_ID = ROIsignals_BN(sub).subj_ID;
    results_BN246(sub).kop = kop;
    results_BN246(sub).mean_kop = MS(sub);
    results_BN246(sub).min_kop = min(kop);
    results_BN246(sub).max_kop = max(kop);
    results_BN246(sub).std_kop = SS(sub);
    results_BN246(sub).cv_kop = CS(sub);
    results_BN246(sub).bins_num = bins_num; 
    results_BN246(sub).entropy_kop = SE(sub);
    results_BN246(sub).sample_failed = sample_failed;
    
end

%% BN210-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_BN');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 210;
    signals = ROIsignals_BN(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = f_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  f_kop2sta(kop, bins_num);
    
    results_BN210(sub).subj_ID = ROIsignals_BN(sub).subj_ID;
    results_BN210(sub).kop = kop;
    results_BN210(sub).mean_kop = MS(sub);
    results_BN210(sub).min_kop = min(kop);
    results_BN210(sub).max_kop = max(kop);
    results_BN210(sub).std_kop = SS(sub);
    results_BN210(sub).cv_kop = CS(sub);
    results_BN210(sub).bins_num = bins_num; 
    results_BN210(sub).entropy_kop = SE(sub);
    results_BN210(sub).sample_failed = sample_failed;
    
end

%% Z1024-ROI signals

% load the ROI signals 
% calculate the MS, SE and Kuramoto order paramater
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_Z');

for sub = 1:295
    sub
    
    % ROI signals
    time_len = 1200;
    node_num = 1024;
    signals = ROIsignals_Z(sub).ROIsignals(1:node_num, 1:time_len);
    
    % z-score normalized signals
    signals_zs = zscore(signals');
    
    %calculate the kuramoto parameter
    kop = f_kop(signals_zs');
    bins_num = 30;
    [MS(sub), SS(sub), CS(sub), SE(sub), sample_failed] =  f_kop2sta(kop, bins_num);
    
    results_Z(sub).subj_ID = ROIsignals_Z(sub).subj_ID;
    results_Z(sub).kop = kop;
    results_Z(sub).mean_kop = MS(sub);
    results_Z(sub).min_kop = min(kop);
    results_Z(sub).max_kop = max(kop);
    results_Z(sub).std_kop = SS(sub);
    results_Z(sub).cv_kop = CS(sub);
    results_Z(sub).bins_num = bins_num; 
    results_Z(sub).entropy_kop = SE(sub);
    results_Z(sub).sample_failed = sample_failed;
    
end

save('step_2_MS_SE_relationship/results.mat','results_HY','results_BN246','results_BN210','results_Z','-v7.3')
