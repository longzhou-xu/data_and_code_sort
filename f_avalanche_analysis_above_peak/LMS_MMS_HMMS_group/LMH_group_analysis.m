%the avalanche analysis in three group: LMS MMS HMS
clc;clear;close all
addpath('CritAnalysisSoftwarePackage2016-04-25');

Smax2fit = 30; 
Smin2fit = 3;
Tmax2fit = 9; 
Tmin2fit=2;

% load the mean synchrony and synchrony entropy of HY96
load(['D:\criticality_cognitive\project_synchrony_avalanche_fl',...
    'uid_Iq\b_MS_SE_STATIC\MS_SE_HY96_RS.mat'])
% select the subject of LMS
[syn_sort,sort_index]=sort(syn);
LMS_Subject(:,1)=syn_sort(1:20,1);
LMS_Subject(:,2)=sort_index(1:20,1);
% avalanche analysis of LMS group
threshold = 1.4;
ava_size_LMS = [];
ava_duration_LMS = []; 
binFirings_LMS = [];
for s=1:length(LMS_Subject(:,2)) % subject
    % load the ROIs' signals
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
        '\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
        num2str(LMS_Subject(s,2)),'.mat']);
    signals_original = rest_HY96_ROI;
    % z-normailiez 
    signal_normalized = z_normalize(signals_original);
    % avalanche analysis (1) : avalanche define
    binsize=1;
    clear avalancheSize binFiring duringTime
    [avalancheSize, binFiring, duringTime, bpSubejct] = avalancheAnalysisBOLD(signal_normalized, threshold, binsize);
    % the LMS group-aggregate size, duration,events in a bin
    ava_size_LMS = [ava_size_LMS; avalancheSize];
    ava_duration_LMS = [ava_duration_LMS; duringTime];
    binFirings_LMS = [binFirings_LMS; binFiring];
end
[bpGroup_LMS]=branchingParameter(binFirings_LMS);
[lifetime_LMS, average_size_LMS]=lifetimeAverageSize(ava_duration_LMS,ava_size_LMS);
for I=5:100
    Kappa_lzGroup_LMS(I,1) = kappa_lz(ava_size_LMS, min(ava_size_LMS), max(ava_size_LMS), I);
end
[alpha_mle_LMS, Smin2fit, Smax2fit] = plmle(ava_size_LMS, 'xmin', Smin2fit, 'xmax', Smax2fit); 
[tau_mle_LMS, Tmin2fit, Tmax2fit] = plmle(ava_duration_LMS,'xmin', Tmin2fit, 'xmax', Tmax2fit);
[gamma_nlsfit_LMS] = createFit_pl_nonlinearLS(lifetime_LMS(Tmin2fit:Tmax2fit), average_size_LMS(Tmin2fit:Tmax2fit));
gamma_fitLMS = gamma_nlsfit_LMS.b;
gamma_mle_LMS = (tau_mle_LMS-1) / (alpha_mle_LMS-1);
save(['avalanches_LMS_',num2str(threshold),'.mat'], ...
    'ava_duration_LMS', 'ava_size_LMS', 'average_size_LMS', 'binFirings_LMS', ...
    'bpGroup_LMS', 'LMS_Subject', 'lifetime_LMS', 'Kappa_lzGroup_LMS', ...
    'alpha_mle_LMS', 'tau_mle_LMS', 'gamma_fitLMS', 'gamma_mle_LMS');

% MMS
[syn_sort,sort_index]=sort(syn);
MMS_Subject(:,1)=syn_sort(201:220,1);
MMS_Subject(:,2)=sort_index(201:220,1);
% avalanche analysis of MMS group
ava_size_MMS = [];
ava_duration_MMS = []; 
binFirings_MMS = [];
for s=1:length(MMS_Subject(:,2)) % subject
    % load the ROIs' signals
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
        '\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
        num2str(MMS_Subject(s,2)),'.mat']);
    signals_original = rest_HY96_ROI;
    % z-normailiez 
    signal_normalized = z_normalize(signals_original);
    % avalanche analysis (1) : avalanche define
    binsize=1;
    clear avalancheSize binFiring duringTime
    [avalancheSize, binFiring, duringTime, bpSubejct] = avalancheAnalysisBOLD(signal_normalized, threshold, binsize);
    % the MMS group-aggregate size, duration,events in a bin
    ava_size_MMS = [ava_size_MMS; avalancheSize];
    ava_duration_MMS = [ava_duration_MMS; duringTime];
    binFirings_MMS = [binFirings_MMS; binFiring];
end
[bpGroup_MMS]=branchingParameter(binFirings_MMS);
[lifetime_MMS, average_size_MMS]=lifetimeAverageSize(ava_duration_MMS,ava_size_MMS);
[alpha_mle_MMS, Smin2fit, Smax2fit] = plmle(ava_size_MMS, 'xmin', Smin2fit, 'xmax', Smax2fit); 
[tau_mle_MMS, Tmin2fit, Tmax2fit] = plmle(ava_duration_MMS,'xmin', Tmin2fit, 'xmax', Tmax2fit);
[gamma_nlsfit_MMS] = createFit_pl_nonlinearLS(lifetime_MMS(Tmin2fit:Tmax2fit), average_size_MMS(Tmin2fit:Tmax2fit));
gamma_fitMMS = gamma_nlsfit_MMS.b;
gamma_mle_MMS = (tau_mle_MMS-1) / (alpha_mle_MMS-1);
for I=5:100
    Kappa_lzGroup_MMS(I,1) = kappa_lz(ava_size_MMS, min(ava_size_MMS), max(ava_size_MMS), I);
end
save(['avalanches_MMS_',num2str(threshold),'.mat'], ...
    "ava_duration_MMS", "ava_size_MMS", "average_size_MMS", "binFirings_MMS", ...
    "bpGroup_MMS", "MMS_Subject", "lifetime_MMS", "Kappa_lzGroup_MMS", ...
    'alpha_mle_MMS', 'tau_mle_MMS', 'gamma_fitMMS', 'gamma_mle_MMS');

% HMS
HMS_Subject(:,1)=syn_sort(276:295,1);
HMS_Subject(:,2)=sort_index(276:295,1);
% avalanche analysis of HMS group
ava_size_HMS = [];
ava_duration_HMS = []; 
binFirings_HMS = [];
for s=1:length(HMS_Subject(:,2)) % subject
    % load the ROIs' signals
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
        '\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
        num2str(HMS_Subject(s,2)),'.mat']);
    signals_original = rest_HY96_ROI;
    % z-normailiez 
    signal_normalized = z_normalize(signals_original);
    % avalanche analysis (1) : avalanche define
    binsize=1;
    clear avalancheSize binFiring duringTime
    [avalancheSize, binFiring, duringTime, bpSubejct] = avalancheAnalysisBOLD(signal_normalized, threshold, binsize);
    % the MMS group-aggregate size, duration,events in a bin
    ava_size_HMS = [ava_size_HMS; avalancheSize];
    ava_duration_HMS = [ava_duration_HMS; duringTime];
    binFirings_HMS = [binFirings_HMS; binFiring];
end
[bpGroup_HMS]=branchingParameter(binFirings_HMS);
[lifetime_HMS, average_size_HMS]=lifetimeAverageSize(ava_duration_HMS,ava_size_HMS);
for I=5:100
    Kappa_lzGroup_HMS(I,1) = kappa_lz(ava_size_HMS, min(ava_size_HMS), max(ava_size_HMS), I);
end
[alpha_mle_HMS, Smin2fit, Smax2fit] = plmle(ava_size_HMS, 'xmin', Smin2fit, 'xmax', Smax2fit); 
[tau_mle_HMS, Tmin2fit, Tmax2fit] = plmle(ava_duration_HMS,'xmin', Tmin2fit, 'xmax', Tmax2fit);
[gamma_nlsfit_HMS] = createFit_pl_nonlinearLS(lifetime_HMS(Tmin2fit:Tmax2fit), average_size_HMS(Tmin2fit:Tmax2fit));
gamma_fitHMS = gamma_nlsfit_HMS.b;
gamma_mle_HMS = (tau_mle_HMS-1) / (alpha_mle_HMS-1);
save(['avalanches_HMS_',num2str(threshold),'.mat'], ...
    "ava_duration_HMS","ava_size_HMS","average_size_HMS","binFirings_HMS",...
    "bpGroup_HMS","HMS_Subject","lifetime_HMS","Kappa_lzGroup_HMS", ...
    'alpha_mle_HMS', 'tau_mle_HMS', 'gamma_fitHMS', 'gamma_mle_HMS');