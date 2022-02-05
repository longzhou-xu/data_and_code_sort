%% step 5
%% power-law analysis of HY96 signals
% group level 
clc;
clear;
close all;
addpath Function
addpath Function/CritAnalysisSoftwarePackage2016-04-25
load(fullfile('step_4_avalanches','results_HY96.mat'),'results_HY96');

for thr = 1:25
    
    sizes = [];
    
    durations = [];
    
    for SUB = 1:295
        
        sizes = [sizes, results_HY96(SUB).avalanches(thr).sta_ava.size];
        durations = [durations, results_HY96(SUB).avalanches(thr).sta_ava.duration];
        
    end
    
    sizes_agg{thr} = sizes;
    
    durations_agg{thr} = durations;
    
end

for thr = 10:2:18 
    % ------------avalanche sizes--------
    thr
    
    sizes = [];
    sizes = sizes_agg{thr};
    
    s_range = [];
    nn = 0;
    x_s = unique(sizes);
    xs_i = find(x_s <= 100 & x_s >= 1);
    for ii = 1 : 10
        for jj = ii + 20 : length(xs_i)
            nn = nn + 1;
            s_range(nn,1) = x_s(xs_i(ii));
            s_range(nn,2) = x_s(xs_i(jj));
        end
    end
    
    % size power law
    powerlaw_group.size(thr).threshold = thr/10;
    powerlaw_group.size(thr).alpha = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_min = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_max = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).p = zeros(length(s_range(:,1)),1);
    
    for num = 1:length(s_range(:,1))
        
        length(s_range(:,1))
        num
        
        % power law fit of size
        [alpha, s_min, s_max] = plmle(sizes, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        % p vaule of power law fit
        p = pvcalc(sizes, alpha, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        
        % save
        powerlaw_group.size(thr).alpha(num) = alpha;
        powerlaw_group.size(thr).s_min(num) = s_min;
        powerlaw_group.size(thr).s_max(num) = s_max;
        powerlaw_group.size(thr).p(num) = p;
        
    end
   
    % ----avalanche durations---

    thr
    
    durations = [];
    durations = durations_agg{thr};
     
    duration_r = [];
    nn = 0;
    x_t = unique(durations);
    xt_i = find(x_t <= 100 & x_t >= 1);
    for ii = 1 : 5
        for jj = ii + 8 : length(xt_i)
            nn = nn + 1;
            duration_r(nn,1) = x_t(xt_i(ii));
            duration_r(nn,2) = x_t(xt_i(jj));
        end
    end 
    
    % duration
    powerlaw_group.duration(thr).threshold = thr/10;
    powerlaw_group.duration(thr).tau = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_min = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_max = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).p = zeros(length(duration_r(:,1)),1);
    
    for num = 1:length(duration_r(:,1))
        
        length(duration_r(:,1))
        num
        
        % power law fit of durations
        [tau, t_min, t_max] = plmle(durations, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        % p value of power law
        p = pvcalc(durations, tau, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        
        % save
        powerlaw_group.duration(thr).tau(num) = tau;
        powerlaw_group.duration(thr).t_min(num) = t_min; 
        powerlaw_group.duration(thr).t_max(num) = t_max;
        powerlaw_group.duration(thr).p(num) = p;
        
    end
    
end
results_HY96(1).powerlaw_group = powerlaw_group;
save(fullfile('step_5_powerlaw_analysis','results_HY96.mat'),'results_HY96','-v7.3'); 

%% power-law analysis of BN246 signals
% group level 
clc;
clear;
close all;
addpath Function
addpath Function/CritAnalysisSoftwarePackage2016-04-25
load(fullfile('step_4_avalanches','results_BN246.mat'),'results_BN246');

for thr = 1:25
    
    sizes = [];
    
    durations = [];
    
    for SUB = 1:295
        
        sizes = [sizes, results_BN246(SUB).avalanches(thr).sta_ava.size];
        durations = [durations, results_BN246(SUB).avalanches(thr).sta_ava.duration];
        
    end
    
    sizes_agg{thr} = sizes;
    
    durations_agg{thr} = durations;
    
end

for thr = 13:2:21 % 13:2:21
    % ----avalanche sizes--------

    thr
    
    sizes = [];
    sizes = sizes_agg{thr};
    
    s_range = [];
    nn = 0;
    x_s = unique(sizes);
    xs_i = find(x_s <= 100 & x_s >= 1);
    for ii = 1 : 10
        for jj = ii + 20 : length(xs_i)
            nn = nn + 1;
            s_range(nn,1) = x_s(xs_i(ii));
            s_range(nn,2) = x_s(xs_i(jj));
        end
    end
    
    % size power law
    powerlaw_group.size(thr).threshold = thr/10;
    powerlaw_group.size(thr).alpha = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_min = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_max = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).p = zeros(length(s_range(:,1)),1);
    
    for num = 1:length(s_range(:,1))
        
        length(s_range(:,1))
        num
        
        % power law fit of size
        [alpha, s_min, s_max] = plmle(sizes, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        % p vaule of power law fit
        p = pvcalc(sizes, alpha, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        
        % save
        powerlaw_group.size(thr).alpha(num) = alpha;
        powerlaw_group.size(thr).s_min(num) = s_min;
        powerlaw_group.size(thr).s_max(num) = s_max;
        powerlaw_group.size(thr).p(num) = p;
        
    end
   
    % ----avalanche durations---

    thr
    
    durations = [];
    durations = durations_agg{thr};
     
    duration_r = [];
    nn = 0;
    x_t = unique(durations);
    xt_i = find(x_t <= 100 & x_t >= 1);
    for ii = 1 : 5
        for jj = ii + 8 : length(xt_i)
            nn = nn + 1;
            duration_r(nn,1) = x_t(xt_i(ii));
            duration_r(nn,2) = x_t(xt_i(jj));
        end
    end 
    
    % duration
    powerlaw_group.duration(thr).threshold = thr/10;
    powerlaw_group.duration(thr).tau = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_min = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_max = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).p = zeros(length(duration_r(:,1)),1);
    
    for num = 1:length(duration_r(:,1))
        
        length(duration_r(:,1))
        num
        
        % power law fit of durations
        [tau, t_min, t_max] = plmle(durations, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        % p value of power law
        p = pvcalc(durations, tau, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        
        % save
        powerlaw_group.duration(thr).tau(num) = tau;
        powerlaw_group.duration(thr).t_min(num) = t_min; 
        powerlaw_group.duration(thr).t_max(num) = t_max;
        powerlaw_group.duration(thr).p(num) = p;
        
    end
    
end
results_BN246(1).powerlaw_group = powerlaw_group;
save(fullfile('step_5_powerlaw_analysis','results_BN246.mat'),'results_BN246','-v7.3'); 

%% power-law analysis of Z1024 signals
% group level 
clc;
clear;
close all;
addpath Function
addpath Function/CritAnalysisSoftwarePackage2016-04-25
load(fullfile('step_4_avalanches','results_Z1024.mat'),'results_Z1024');

sub_list = [1:5,7:295];
for thr = 1:25
    
    sizes = [];
    
    durations = [];
    
    for S = 1:length(sub_list)
        
        SUB = sub_list(S);
        
        sizes = [sizes, results_Z1024(SUB).avalanches(thr).sta_ava.size];
        durations = [durations, results_Z1024(SUB).avalanches(thr).sta_ava.duration];
        
    end
    
    sizes_agg{thr} = sizes;
    
    durations_agg{thr} = durations;
    
end

% ------------avalanche sizes--------

for thr = 16:2:24
    thr
    
    sizes = [];
    sizes = sizes_agg{thr};
    
    s_range = [];
    nn = 0;
    x_s = unique(sizes);
    xs_i = find(x_s <= 100 & x_s >= 1);
    for ii = 1 : 10
        for jj = ii + 20 : length(xs_i)
            nn = nn + 1;
            s_range(nn,1) = x_s(xs_i(ii));
            s_range(nn,2) = x_s(xs_i(jj));
        end
    end
    
    % size power law
    powerlaw_group.size(thr).threshold = thr/10;
    powerlaw_group.size(thr).alpha = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_min = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).s_max = zeros(length(s_range(:,1)),1);
    powerlaw_group.size(thr).p = zeros(length(s_range(:,1)),1);
    
    for num = 1:length(s_range(:,1))
        
        length(s_range(:,1))
        num
        
        % power law fit of size
        [alpha, s_min, s_max] = plmle(sizes, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        % p vaule of power law fit
        p = pvcalc(sizes, alpha, 'xmin', s_range(num,1), 'xmax', s_range(num,2));
        
        % save
        powerlaw_group.size(thr).alpha(num) = alpha;
        powerlaw_group.size(thr).s_min(num) = s_min;
        powerlaw_group.size(thr).s_max(num) = s_max;
        powerlaw_group.size(thr).p(num) = p;
        
    end
   
% ----avalanche durations---

    thr
    
    durations = [];
    durations = durations_agg{thr};
     
    duration_r = [];
    nn = 0;
    x_t = unique(durations);
    xt_i = find(x_t <= 100 & x_t >= 1);
    for ii = 1 : 5
        for jj = ii + 8 : length(xt_i)
            nn = nn + 1;
            duration_r(nn,1) = x_t(xt_i(ii));
            duration_r(nn,2) = x_t(xt_i(jj));
        end
    end 
    
    % duration
    powerlaw_group.duration(thr).threshold = thr/10;
    powerlaw_group.duration(thr).tau = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_min = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).t_max = zeros(length(duration_r(:,1)),1);
    powerlaw_group.duration(thr).p = zeros(length(duration_r(:,1)),1);
    
    for num = 1:length(duration_r(:,1))
        
        length(duration_r(:,1))
        num
        
        % power law fit of durations
        [tau, t_min, t_max] = plmle(durations, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        % p value of power law
        p = pvcalc(durations, tau, 'xmin', duration_r(num,1), 'xmax', duration_r(num,2));
        
        % save
        powerlaw_group.duration(thr).tau(num) = tau;
        powerlaw_group.duration(thr).t_min(num) = t_min; 
        powerlaw_group.duration(thr).t_max(num) = t_max;
        powerlaw_group.duration(thr).p(num) = p;
        
    end
    
end
results_Z1024(1).powerlaw_group = powerlaw_group;
save(fullfile('step_5_powerlaw_analysis','results_Z1024.mat'),'results_Z1024','-v7.3'); 