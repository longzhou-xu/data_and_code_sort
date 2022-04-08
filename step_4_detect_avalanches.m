%% step 4: define the avalanche
clc;
clear;
close all
addpath('Function')

%% detect the avalanches of HY96-ROIsignals
clc;
clear;
close all;
timebinsize = 1;
ROI_number = 96;
time_length = 1200;
for SUB = 1 : 295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    load(fullfile('step_3_events', 'ROI_level', 'HY96', [subj_ID, '.mat']), 'peakevents');
    for THR = 1: 25
        raster = peakevents.ithr(THR).raster;
        avalanches_HY96(SUB).threshold(THR).threshold = THR/10;
        avalanches_HY96(SUB).threshold(THR).sta_ava = xlz_avalanches(raster, ROI_number, time_length, timebinsize); 
    end
end
% save
save(fullfile('step_4_avalanches','avalanches_HY96.mat'), 'avalanches_HY96','-v7.3');

%% detect the avalanches of BN246-ROIsignals
clc;
clear;
close all;
timebinsize = 1;
ROI_number = 246;
time_length = 1200;
for SUB = 1 : 295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    load(fullfile('step_3_events', 'ROI_level', 'BN246', [subj_ID, '.mat']), 'peakevents');
    for THR = 1: 25
         raster = peakevents.ithr(THR).raster;
        avalanches_BN246(SUB).threshold(THR).threshold = THR/10;
        avalanches_BN246(SUB).threshold(THR).sta_ava = xlz_avalanches(raster, ROI_number, time_length, timebinsize); 
    end
end
% save
save(fullfile('step_4_avalanches','avalanches_BN246.mat'), 'avalanches_BN246','-v7.3');

%% detect the avalanches of BN210-ROIsignals
clc;
clear;
close all;
timebinsize = 1;
ROI_number = 210;
time_length = 1200;
for SUB = 1 : 295
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    load(fullfile('step_3_events', 'ROI_level', 'BN210', [subj_ID, '.mat']), 'peakevents');
    for THR = 1: 25
            raster = peakevents.ithr(THR).raster;
            avalanches_BN210(SUB).threshold(THR).threshold = THR/10;
            avalanches_BN210(SUB).threshold(THR).sta_ava = xlz_avalanches(raster, ROI_number, time_length, timebinsize); 
    end
end
% save
save(fullfile('step_4_avalanches', 'avalanches_BN210.mat'), 'avalanches_BN210', '-v7.3');

%% detect the avalanches of Z1024-ROIsignals
clc;
clear;
close all;
timebinsize = 1;
ROI_number = 1024;
time_length = 1200;
sub_list = [1:5,7:295];
for S = 1 : length(sub_list)
    SUB = sub_list(S);
    subj_ID = ['sub_', num2str(SUB, '%03d')]
    load(fullfile('step_3_events', 'ROI_level', 'Z1024', [subj_ID, '.mat']), 'peakevents');
    for THR = 1: 25
        raster = peakevents.ithr(THR).raster;
        avalanches_Z1024(SUB).threshold(THR).threshold = THR/10;
        avalanches_Z1024(SUB).threshold(THR).sta_ava = xlz_avalanches(raster, ROI_number, time_length, timebinsize); 
    end
end
% save
save(fullfile('step_4_avalanches', 'avalanches_Z1024.mat'), 'avalanches_Z1024', '-v7.3');
