%% step three: show the MS-SE relationship
clc;
clear;
close all
%% setting environment
addpath('Function')
%% setting data path:
subj_folder = 'step_3_events/ROI_level/HY96/';
%% detect the avalanches
timebinsize = 1;
time_length = 1200;
for SUB = 1 : 295
    file_path = [subj_folder, 'sub',num2str(SUB,'%.3d'),'_peakevents.mat'];
    load(file_path);
    for thr = 1: 25
        raster = [];
        raster = events.ithr(thr).raster;
        ava.isub(SUB).ithr(thr) = detect_avalanches(raster, time_length, timebinsize); 
    end
end
load('step_1_signals/subject.mat', 'subject');
ava.subjectID = subject;
ava.threshold = 0.1 : 0.1 : 2.5;
%% saving
mkdir('step_4_avalanche/ROI_level/HY96');
save('step_4_avalanche/ROI_level/HY96/peakevents_aval.mat','ava');

