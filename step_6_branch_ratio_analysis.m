%% step three: show the MS-SE relationship
clc;
clear;
close all
%% setting environment
addpath('Function');
addpath('Function\cifti-matlab-master');
%% setting data path:
subj_folder = 'step_4';
subj_dir = dir([subj_folder,'\sub*']);
%% collect the set of size and duration 
for SUB = 1:length(subj_dir)
    file_path = ['step_4\',subj_dir(SUB).name,'\avalanches.mat']
    load(file_path);
    for thr = 1 : 40 % 0.1:0.1:4.0
        for Run = 1 : 4
            % 
            shape = ava.ithr(thr).irun(Run).ava.shape;
        end
    end
end