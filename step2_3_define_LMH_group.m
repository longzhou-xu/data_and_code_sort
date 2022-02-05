%% define the LMS、MMS and HMS groups
clc;
clear;
close all;

load(fullfile('step_2_MS_SE_relationship','results.mat'));

for sub = 1:295
    % import the mean synchronization
    MS(sub) = results_HY(sub).mean_kop;
    subj_ID(sub) = results_HY(sub).subj_ID;
end

% sort the MS in ascending order
[MS_value, MS_index] = sort(MS);

% LMS group
LMHgroup.LMS(:,1) = subj_ID(MS_index(1:20));
LMHgroup.LMS(:,2) = MS_value(1:20);
LMHgroup.LMS(:,3) = MS_index(1:20);

% MMS group
LMHgroup.MMS(:,1) = subj_ID(MS_index(200:219));
LMHgroup.MMS(:,2) = MS_value(200:219);
LMHgroup.MMS(:,3) = MS_index(200:219);

% HMS group
LMHgroup.HMS(:,1) = subj_ID(MS_index(276:295));
LMHgroup.HMS(:,2) = MS_value(276:295);
LMHgroup.HMS(:,3) = MS_index(276:295);

save(fullfile('step_2_MS_SE_relationship','results.mat'),'results_HY','results_BN210','results_BN246','results_Z','LMHgroup');