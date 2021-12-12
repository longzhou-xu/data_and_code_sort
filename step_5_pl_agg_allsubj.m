%% step 5 power-law analysis at all-subject-aggragate level
%% step 5 - 1 aggreagte the size and duration
clc;
clear;
close all;
% set
addpath Function
addpath Function/CritAnalysisSoftwarePackage2016-04-25
load('step_4_avalanche/ROI_level/HY96/peakevents_aval.mat','ava');
for thr = 1:25
    sizes = [];
    durations = [];
    for SUB = 1:295
        sizes = [sizes,ava.isub(SUB).ithr(thr).size];
        durations = [durations,ava.isub(SUB).ithr(thr).duration];
    end

    mkdir('step_5_pl_analy/ROI_level/HY96/agg_l_st');
    save(['step_5_pl_analy/ROI_level/HY96/agg_l_st/st_',num2str(thr),'.mat'], ...
        'sizes','durations');
end
%% step 5 - 2 power-law analysis
clc;
clear;
close all;
addpath Function
addpath Function/CritAnalysisSoftwarePackage2016-04-25
for thr = 1:25
    thr
    sizes = [];
    durations = [];
    load(['step_5_pl_analy/ROI_level/HY96/agg_l_st/st_',num2str(thr),'.mat']);
    % ------------------avalanche size----
    size_range = [];
    nn = 0;
    x_s = unique(sizes);
    xs_i = find(x_s <= 100 & x_s >= 1);
    for ii = 1 : length(xs_i)-10
        for jj = ii + 30 : length(xs_i)
            nn = nn + 1;
            size_range(nn,1) = x_s(xs_i(ii));
            size_range(nn,2) = x_s(xs_i(jj));
        end
    end
    % size power law
    spl.alpha = zeros(length(size_range(:,1)),1);
    spl.s_min = zeros(length(size_range(:,1)),1);
    spl.s_max = zeros(length(size_range(:,1)),1);
    spl.p = zeros(length(size_range(:,1)),1);
    for num = 1:length(size_range(:,1))
        num
        tic;
        [alpha, s_min, s_max] = ...
            plmle(sizes, 'xmin', size_range(num,1), 'xmax', size_range(num,2));
        p = pvcalc(sizes, spl.alpha(num), ...
            'xmin', size_range(num,1), 'xmax', size_range(num,2));
        spl.alpha(num) = alpha;
        spl.s_min(num) = s_min;
        spl.s_max(num) = s_max;
        spl.p(num) = p;
        toc;
    end
    % ----avalanche duration---
    duration_range = [];
    nn = 0;
    x_t = unique(durations);
    xt_i = find(x_t <= 100 & x_t >= 1);
    for ii = 1 : length(xt_i)-10
        for jj = ii + 10 : length(xt_i)
            nn = nn + 1;
            duration_range(nn,1) = x_t(xt_i(ii));
            duration_range(nn,2) = x_t(xt_i(jj));
        end
    end 
    % duration
    tpl.tau = zeros(length(duration_range(:,1)),1);
    tpl.t_min = zeros(length(duration_range(:,1)),1);
    tpl.t_max = zeros(length(duration_range(:,1)),1);
    tpl.p = zeros(length(duration_range(:,1)),1);
    for num = 1:length(duration_range(:,1))
        [tau, t_min, t_max] = plmle(durations, 'xmin', ...
            duration_range(num,1), 'xmax', duration_range(num,2));
        p = pvcalc(durations, tpl.tau(num), ...
            'xmin', duration_range(num,1), 'xmax', duration_range(num,2));
        tpl.tau(num) = tau;
        tpl.t_min(num) = t_min; 
        tpl.t_max(num) = t_max;
        tpl.p(num) = p;
    end
    mkdir('step_5_pl_analy/ROI_level/HY96/agg_pl');
    save(['step_5_pl_analy/ROI_level/HY96/agg_pl/pl_',num2str(thr),'.mat'], ...
        'spl','tpl');
end

    
    