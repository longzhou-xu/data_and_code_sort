%% step three: detect the peak events
clc;
clear;
close all
addpath('Function');
%% detect the peak events of HY-96 ROI signals
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'), 'ROIsignals_HY');
mkdir(fullfile('step_3_events', 'ROI_level', 'HY96'));
for SUB = 1:295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    II = 0;
    clear peakevents
    for threshold = 0.1 : 0.1 : 3
        II = II + 1;
        signals = ROIsignals_HY(SUB).ROIsignals;
        signals = zscore(signals')';
        time_length = 1200;
        node_num = 96;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = xlz_peakevents(signals, threshold, time_length, node_num);
        peakevents.ithr(II).threshold = threshold;
    end
    peakevents.threshold = 0.1 : 0.1 : 3;
    % save
    save(fullfile('step_3_events', 'ROI_level','HY96', ['sub_', num2str(SUB,'%03d'),'.mat']), 'peakevents');
end
%% detect the peak events of BN-246 ROI signals
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'), 'ROIsignals_BN');
mkdir(fullfile('step_3_events','ROI_level','BN246'));
for SUB = 1:295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    II = 0;
    clear peakevents
    for threshold = 0.1 : 0.1 : 3
        II = II + 1;
        signals = ROIsignals_BN(SUB).ROIsignals;
        signals = zscore(signals')';
        time_length = 1200;
        node_num = 246;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = xlz_peakevents(signals, threshold, time_length, node_num);
        peakevents.ithr(II).threshold = threshold;
    end
    peakevents.threshold = 0.1 : 0.1 : 3;
    % save
    save(fullfile('step_3_events','ROI_level','BN246',['sub_', num2str(SUB,'%03d'),'.mat']), 'peakevents');
end

%% detect the peak events of BN-210 ROI signals
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'), 'ROIsignals_BN');
mkdir(fullfile('step_3_events','ROI_level','BN210'));
for SUB = 1:295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    II = 0;
    clear peakevents
    for threshold = 0.1 : 0.1 : 3
        II = II + 1;
        signals = ROIsignals_BN(SUB).ROIsignals;
        signals = zscore(signals')';
        signals = signals(1:210, :);
        time_length = 1200;
        node_num = 210;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = xlz_peakevents(signals, threshold, time_length, node_num);
        peakevents.ithr(II).threshold = threshold;
    end
    peakevents.threshold = 0.1 : 0.1 : 3;
    % save
    save(fullfile('step_3_events','ROI_level','BN210',['sub_', num2str(SUB,'%03d'),'.mat']), 'peakevents');
end

%% detect the peak events of Z-1024 ROI signals
load(fullfile('step_1_extract_ROI_signals', 'ROI_signals', 'ROIsignals.mat'), 'ROIsignals_Z');
mkdir(fullfile('step_3_events', 'ROI_level', 'Z1024'));
sub_list = [1:5,7:295];
for s = 1:length(sub_list)
    SUB = sub_list(s);
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    II = 0;
    clear peakeventsZ
    for threshold = 0.1 : 0.1 : 3
        II = II + 1;
        signals = ROIsignals_Z(SUB).ROIsignals;
        signals = zscore(signals')';
        time_length = 1200;
        node_num = 1024;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = xlz_peakevents(signals, threshold, time_length, node_num);
    end
    peakevents.threshold = 0.1 : 0.1 : 3;
    % save
    save(fullfile('step_3_events', 'ROI_level', 'Z1024', ['sub_', num2str(SUB, '%03d'),'.mat']), 'peakevents');
end
