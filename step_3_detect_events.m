%% step three: show the MS-SE relationship
clc;
clear;
close all
addpath('Function');

%% detect the peak events of HY-96 ROI signals
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_HY');
load(fullfile('step_2_MS_SE_relationship','results.mat'), 'results_HY');
mkdir(fullfile('step_3_events','ROI_level','HY96'));
for sub = 1:295
    sub
    II = 0;
    clear peakevents peakEvents
    for threshold = 0.1 : 0.1 : 3
        
        II = II + 1;
        signals = ROIsignals_HY(sub).ROIsignals;
        signals = zscore(signals')';
        
        time_length = 1200;
        node_num = 96;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = ...
                f_peakevents(signals, threshold, time_length, node_num);
            
        peakEvents(II).threshold = threshold;
        peakEvents(II).fingerprint = peakevents.ithr(II).fingerprint;
        
    end
    
    peakevents.threshold = 0.1 : 0.1 : 3;
    
    save(fullfile('step_3_events','ROI_level','HY96',['sub_', num2str(sub,'%03d'),'.mat']), 'peakevents');
    
    results_HY(sub).peakEvents = peakEvents;
end

%% detect the peak events of BN-246 ROI signals
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_BN');
load(fullfile('step_2_MS_SE_relationship','results.mat'), 'results_BN246');
mkdir(fullfile('step_3_events','ROI_level','BN246'));
for sub = 1:295
    sub
    II = 0;
    clear peakevents peakEvents
    for threshold = 0.1 : 0.1 : 3
        
        II = II + 1;
        signals = ROIsignals_BN(sub).ROIsignals;
        signals = zscore(signals')';
        
        time_length = 1200;
        node_num = 246;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = ...
                f_peakevents(signals, threshold, time_length, node_num);
            
        peakEvents(II).threshold = threshold;
        peakEvents(II).fingerprint = peakevents.ithr(II).fingerprint;
    end
    
    peakevents.threshold = 0.1 : 0.1 : 3;
    
    save(fullfile('step_3_events','ROI_level','BN246',['sub_', num2str(sub,'%03d'),'.mat']), 'peakevents');
    
    results_BN246(sub).peakEvents = peakEvents;
end

%% detect the peak events of BN-210 ROI signals
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_BN');
load(fullfile('step_2_MS_SE_relationship','results.mat'), 'results_BN210');
mkdir(fullfile('step_3_events','ROI_level','BN210'));
for sub = 1:295
    sub
    II = 0;
    clear peakevents peakEvents
    for threshold = 0.1 : 0.1 : 3
        
        II = II + 1;
        signals = ROIsignals_BN(sub).ROIsignals;
        signals = zscore(signals')';
        signals = signals(1:210, :);
        
        time_length = 1200;
        node_num = 210;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = ...
                f_peakevents(signals, threshold, time_length, node_num);
            
        peakEvents(II).threshold = threshold;
        peakEvents(II).fingerprint = peakevents.ithr(II).fingerprint;
    end
    
    peakevents.threshold = 0.1 : 0.1 : 3;
    
    save(fullfile('step_3_events','ROI_level','BN210',['sub_', num2str(sub,'%03d'),'.mat']), 'peakevents');
    
    results_BN210(sub).peakEvents = peakEvents;
end

%% detect the peak events of Z-1024 ROI signals
load(fullfile('step1_signals','ROI_signals','ROIsignals.mat'), 'ROIsignals_Z');
load(fullfile('step_2_MS_SE_relationship','results.mat'), 'results_Z');
mkdir(fullfile('step_3_events','ROI_level','Z1024'));
sub_list = [1:5,7:295];
for s = 1:length(sub_list)
    sub = sub_list(s);
    sub
    II = 0;
    clear peakeventsZ peakEvents
    for threshold = 0.1 : 0.1 : 3
        
        II = II + 1;
        signals = ROIsignals_Z(sub).ROIsignals;
        signals = zscore(signals')';
        
        time_length = 1200;
        node_num = 1024;
        [peakevents.ithr(II).raster, peakevents.ithr(II).fingerprint] = ...
                f_peakevents(signals, threshold, time_length, node_num);
            
        peakEvents(II).threshold = threshold;
        peakEvents(II).fingerprint = peakevents.ithr(II).fingerprint;
    end
    
    peakevents.threshold = 0.1 : 0.1 : 3;
    
    save(fullfile('step_3_events','ROI_level','Z1024',['sub_', num2str(sub,'%03d'),'.mat']), 'peakevents');
    
    results_Z(sub).peakEvents = peakEvents;
end
%% save
load(fullfile('step_2_MS_SE_relationship','results.mat'), 'LMHgroup');
save(fullfile('step_3_events','results.mat'), 'results_HY', 'results_BN246', 'results_BN210', 'results_Z', 'LMHgroup','-v7.3');