%% step three: show the MS-SE relationship
clc;
clear;
close all
addpath('Function')

%% detect the avalanches of HY96-ROIsignals
load(fullfile('step_3_events','results.mat'), 'results_HY');
timebinsize = 1;
channel_length = 96;
time_length = 1200;
for SUB = 1 : 295
    SUB
    
    for thr = 1: 25
        fingerprint = results_HY(SUB).peakEvents(thr).fingerprint;
        results_HY(SUB).avalanches(thr).threshold = thr/10;
        results_HY(SUB).avalanches(thr).sta_ava = f_avalanches(fingerprint, channel_length, time_length, timebinsize); 
    end
    
end
% save
results_HY96 = results_HY;
load(fullfile('step_3_events','results.mat'), 'LMHgroup');
results_HY96(1).LMHgroups = LMHgroup;
save(fullfile('step_4_avalanches','results_HY96.mat'),'results_HY96','-v7.3');

%% detect the avalanches of BN246-ROIsignals
load(fullfile('step_3_events','results.mat'), 'results_BN246');
timebinsize = 1;
channel_length = 246;
time_length = 1200;
for SUB = 1 : 295
    SUB
    
    for thr = 1: 25
        fingerprint = results_BN246(SUB).peakEvents(thr).fingerprint;
        results_BN246(SUB).avalanches(thr).threshold = thr/10;
        results_BN246(SUB).avalanches(thr).sta_ava = f_avalanches(fingerprint, channel_length, time_length, timebinsize); 
    end
    
end
% save
save(fullfile('step_4_avalanches','results_BN246.mat'),'results_BN246','-v7.3');

%% detect the avalanches of BN210-ROIsignals
load(fullfile('step_3_events','results.mat'), 'results_BN210');
timebinsize = 1;
channel_length = 210;
time_length = 1200;
for SUB = 1 : 295
    SUB
    
    for thr = 1: 25
        fingerprint = results_BN210(SUB).peakEvents(thr).fingerprint;
        results_BN210(SUB).avalanches(thr).threshold = thr/10;
        results_BN210(SUB).avalanches(thr).sta_ava = f_avalanches(fingerprint, channel_length, time_length, timebinsize); 
    end
    
end
% save
save(fullfile('step_4_avalanches','results_BN210.mat'),'results_BN210','-v7.3');

%% detect the avalanches of Z1024-ROIsignals
load(fullfile('step_3_events','results.mat'), 'results_Z');
timebinsize = 1;
channel_length = 1024;
time_length = 1200;
sub_list = [1:5,7:295];
for S = 1 : length(sub_list)
    
    SUB = sub_list(S)
    
    for thr = 1: 25
        fingerprint = results_Z(SUB).peakEvents(thr).fingerprint;
        results_Z(SUB).avalanches(thr).threshold = thr/10;
        results_Z(SUB).avalanches(thr).sta_ava = f_avalanches(fingerprint, channel_length, time_length, timebinsize); 
    end
    
end
% save
results_Z1024 = results_Z;
save(fullfile('step_4_avalanches','results_Z1024.mat'),'results_Z1024','-v7.3');