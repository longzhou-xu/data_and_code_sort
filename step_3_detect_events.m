%% step three: show the MS-SE relationship
clc;
clear;
close all
%% setting environment
addpath('Function');
%% setting data path:
signals_folder = 'step_1_signals/ROI_signals/HY_96';
mkdir('step_3_events/ROI_level/HY96');
%% detect the event
for SUB = 1:295
    file_path = [signals_folder,'/sub',num2str(SUB),'.mat'];
    load(file_path);
    II = 0;
    for threshold = 0.1 : 0.1 : 3
        II = II + 1;
        signals = rest_HY96_ROI;
        signals = zscore(signals)';
        time_length = 1200;
        node_num = 96;
        [events.ithr(II).raster, events.ithr(II).fingerprint] = ...
                detect_events(signals, threshold, time_length, node_num);
    end
    events.threshold = 0.1 : 0.1 : 3;
    file_path = ['step_3_events/ROI_level/HY96/sub',num2str(SUB,'%.3d'),'_peakevents.mat'];
    save(file_path,'events');
end