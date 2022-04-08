%% convet the avalanche to branching process

%% HY96
clc; clear; close all;
addpath('Function');
load(fullfile('step_4_avalanches', 'avalanches_HY96.mat'))
sub_list = 1:295;
for N = 1:length(sub_list)
    sub = sub_list(N)
    threshold = 0;
    for TT = 1:25
        threshold = threshold + 0.1;
        shape = avalanches_HY96(sub).threshold(TT).sta_ava.shape;
        BP_HY96.sub(sub).thr(TT).branching_ratio_series = xlz_shape2branchingprocess(shape);
        BP_HY96.sub(sub).thr(TT).threshold = threshold;
    end
end
mkdir('step_6_branching_process_analysis');
save(fullfile('step_6_branching_process_analysis', 'branching_process_HY96.mat'), 'BP_HY96');

%% BN246
clc; clear; close all;
addpath('Function');
load(fullfile('step_4_avalanches', 'avalanches_BN246.mat'))
sub_list = 1:295;
for N = 1:length(sub_list)
    sub = sub_list(N)
    threshold = 0;
    for TT = 1:25
        threshold = threshold + 0.1;
        shape = avalanches_BN246(sub).threshold(TT).sta_ava.shape;
        BP_BN246.sub(sub).thr(TT).branching_ratio_series = xlz_shape2branchingprocess(shape);
        BP_BN246.sub(sub).thr(TT).threshold = threshold;
    end
end
mkdir('step_6_branching_process_analysis');
save(fullfile('step_6_branching_process_analysis', 'branching_process_BN246.mat'), 'BP_BN246');

%% BN210
clc; clear; close all;
addpath('Function');
load(fullfile('step_4_avalanches', 'avalanches_BN210.mat'))
sub_list = 1:295;
for N = 1:length(sub_list)
    sub = sub_list(N)
    threshold = 0;
    for TT = 1:25
        threshold = threshold + 0.1;
        shape = avalanches_BN210(sub).threshold(TT).sta_ava.shape;
        BP_BN210.sub(sub).thr(TT).branching_ratio_series = xlz_shape2branchingprocess(shape);
        BP_BN210.sub(sub).thr(TT).threshold = threshold;
    end
end
mkdir('step_6_branching_process_analysis');
save(fullfile('step_6_branching_process_analysis', 'branching_process_BN210.mat'), 'BP_BN210');

%% Z1024
clc; clear; close all;
addpath('Function');
load(fullfile('step_4_avalanches', 'avalanches_Z1024.mat'))
sub_list = [1:5, 7:295];
for N = 1:length(sub_list)
    sub = sub_list(N)
    threshold = 0;
    for TT = 1:25
        threshold = threshold + 0.1;
        shape = avalanches_Z1024(sub).threshold(TT).sta_ava.shape;
        BP_Z1024.sub(sub).thr(TT).branching_ratio_series = xlz_shape2branchingprocess(shape);
        BP_Z1024.sub(sub).thr(TT).threshold = threshold;
    end
end
mkdir('step_6_branching_process_analysis');
save(fullfile('step_6_branching_process_analysis', 'branching_process_Z1024.mat'), 'BP_Z1024');
