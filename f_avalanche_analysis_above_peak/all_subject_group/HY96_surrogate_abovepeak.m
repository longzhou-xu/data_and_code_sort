%% The first part: avalanche define
clc;clear;close all
addpath('CritAnalysisSoftwarePackage2016-04-25');
I_T = 0;
THD=1:0.1:2.6;
for Threshold = 1 : 0.1 : 2.6 
    I_T = I_T + 1;
    AvaS_AllSubject = [];
    AvaD_AllSubject = []; 
    BinE_AllSubject = [];
    for S=1:295 % subject
        disp(['subject=',num2str(S)]);
       %% load the ROIs' signals and normalize the signals
        load(['D:\supplyment_dataANDcode\l_surrogate_data\power_fit_test\HY96\surrogate_data_noremain_FC\surrogate_sub_',num2str(S),'.mat']); %load the ROI signals
        Signals_original = Surrogate_data;
        % z-normailiez 
        Signal_normalized = z_normalize(Signals_original);
       %% Avalanche analysis
        % avalanche define and ststistic
        BinSize=1;
        [Event{S,1}, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
        % the group-aggregate size, duration, events
        AvaS_AllSubject = [AvaS_AllSubject; AvaS{S,I_T}];
        AvaD_AllSubject = [AvaD_AllSubject; DurT{S,I_T}];
        BinE_AllSubject = [BinE_AllSubject; BinE{S,I_T}];
    end
    save(['Avalanches_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_AllSubject', 'AvaD_AllSubject', 'BinE_AllSubject', 'Event');
end
save('Avalanches_define.mat', 'BinE', 'AvaS', 'DurT', 'AvaSH', 'AvaDH', 'BP', 'BP_other');
%% power law analysis
clc; clear; close all
load('Avalanches_1.4.mat')
SizeFitRange = [6,30];
DurationFitRange = [4,8];
[Avalanche_surrogate_nofc, Criterion_surrogate_nofc, BP_GROUP_surrogate_nofc] = ...
    CriticalAnalysis_powerlaw(AvaD_AllSubject, AvaS_AllSubject, ...
    BinE_AllSubject, SizeFitRange, DurationFitRange);
% plot the bp
load('D:\supplyment_dataANDcode\c_MS_SE_static\MS_SE_HY96_RS.mat', 'syn')
load('Avalanches_define.mat', 'BP_other')
figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name',  'Power law analysis', ...
      'Position', [0.1 0.3 0.5 0.5]);
AX1 = subplot(1,1,1);
plot(syn, BP_other(:,1,1) ,...
    'Color', [0.50, 0.50, 0.50], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
[R,P] = corr(syn, BP_other(:,1,1));
xlabel('MS');
ylabel('BP')
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
save('surrogate_data_noremain_FC_avalanche_analysis.mat')
%% power law analysis robustness
clc; clear; close all 
addpath('D:\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25')
load('Avalanches_1.4.mat')
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_AllSubject, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_AllSubject, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end

[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_AllSubject, AvaS_AllSubject);
for Dura_min = 1 : 5
    for Dura_max = 9 : 20
        X = Dura_min
        Y = Dura_max-8
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_AllSubject,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_AllSubject, Tau(X, Y), 'xmin',Tmin2fit(X, Y),...
            'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), ...
            AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('HY96_surrogate_powerlaw_diff_size_duration_1.6_allsubject.mat');
%%
clc;clear;close all
load('HY96_surrogate_powerlaw_diff_size_duration_1.6_allsubject.mat')
Alpha_r = reshape(Alpha,I*J,1);
P_alpha_r = reshape(P_alpha,I*J,1);
Smin2fit_r = reshape(Smin2fit,I*J,1);
Smax2fit_r = reshape(Smax2fit,I*J,1);
Alpha_rc = Alpha_r(find(P_alpha_r>0.2));
[Alpha_rc,Alpha_rc_index] = sort(Alpha_rc);
P_alpha_rc = P_alpha_r(find(P_alpha_r>0.2));
P_alpha_rc = P_alpha_rc(Alpha_rc_index);
Smin2fit_rc = Smin2fit_r(find(P_alpha_r>0.2));
Smin2fit_rc = Smin2fit_rc(Alpha_rc_index);
Smax2fit_rc =Smax2fit_r(find(P_alpha_r>0.2));
Smax2fit_rc = Smax2fit_rc(Alpha_rc_index);

Tau_r = reshape(Tau,X*Y,1);
P_tau_r = reshape(P_tau,X*Y,1);
Tmin2fit_r = reshape(Tmin2fit,X*Y,1);
Tmax2fit_r = reshape(Tmax2fit,X*Y,1);
Gamma_fit_r = reshape(Gamma_fit,X*Y,1);
Gamma_fit_a_r = reshape(Gamma_fit_a,X*Y,1);
Tau_rc = Tau_r(find(P_tau_r>0.2));
[Tau_rc,Tau_rc_index] = sort(Tau_rc);
P_tau_rc = P_tau_r(find(P_tau_r>0.2));
P_tau_rc = P_tau_rc(Tau_rc_index);
Tmin2fit_rc = Tmin2fit_r(find(P_tau_r>0.2));
Tmin2fit_rc = Tmin2fit_rc(Tau_rc_index);
Tmax2fit_rc = Tmax2fit_r(find(P_tau_r>0.2));
Tmax2fit_rc = Tmax2fit_rc(Tau_rc_index);
Gamma_fit_rc = Gamma_fit_r(find(P_tau_r>0.2));
Gamma_fit_rc = Gamma_fit_rc(Tau_rc_index);
Gamma_fit_a_rc = Gamma_fit_a_r(find(P_tau_r>0.2));
Gamma_fit_a_rc = Gamma_fit_a_rc(Tau_rc_index);

for M = 1 : length(Alpha_rc)
    for N = 1 : length(Tau_rc) 
        Gamma_alpha_tau(M,N) = (Tau_rc(N)-1) / (Alpha_rc(M)-1);
        Gamma_fit_1(M,N) = Gamma_fit_rc(N);
    end
end
DLAT = abs(Gamma_alpha_tau-Gamma_fit_1);

AX1 = subplot(3,2,1);
imagesc([30,60],[1,10],Alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\alpha$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,2);
imagesc([30,60],[1,10],P_alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,4);
imagesc([9 20],[1 4],Tau);
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\tau$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,5);
imagesc([9 20],[1 4],P_tau)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,6);
imagesc([9 20],[1 4],Gamma_fit)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\gamma$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');

AX1 = subplot(3,2,5);
Gamma_alpha_tau_1 = Gamma_alpha_tau;
for X = 1 : 107
    for Y = 1 : 1
        if Gamma_alpha_tau(X,Y) >= max(max(Gamma_fit))+2
           Gamma_alpha_tau_1(X,Y) = max(max(Gamma_fit))+2;
        end
    end
end
imagesc(Gamma_alpha_tau_1)
% xticks([2 4 6 8])
% xticklabels({num2str(Tau_rc(2)),num2str(Tau_rc(4)),num2str(Tau_rc(6)),...
%     num2str(Tau_rc(8))});
% yticks([1 50 100 150 200])
% yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
%     num2str(Alpha_rc(150)),num2str(Alpha_rc(200))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\frac{1-\tau}{1-\alpha}$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,6);
for X = 1 : 107
    for Y = 1 : 1
        if DLAT(X,Y) <= 0.1
           Delta(X,Y) = DLAT(X,Y);
        end
        if DLAT(X,Y) > 0.1
            Delta(X,Y) = 0.15;
        end
    end
end
imagesc(Delta)
% xticks([2 4 6 8])
% xticklabels({num2str(Tau_rc(2)),num2str(Tau_rc(4)),num2str(Tau_rc(6)),...
%     num2str(Tau_rc(8))});
% yticks([1 50 100 150 200])
% yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
%     num2str(Alpha_rc(150)),num2str(Alpha_rc(200))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$|\frac{1-\tau}{1-\alpha}-\gamma|$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
%% plot the event point figure of HY-96 surrogate

clc;clear;close all
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\HY_96_2');
load('Avalanches_1.4.mat', 'Event');% load the data of event
load('D:\supplyment_dataANDcode\c_MS_SE_static\MS_SE_HY96_RS.mat', 'syn');%load the data of synchronize
[MS_sorted, Index_MS_sorted] = sort(syn); % rank the mean synchronize
LMS(:,1) = MS_sorted(1:20,1);
LMS(:,2) = Index_MS_sorted(1:20,1);
MMS(:,1) = MS_sorted(201:220,1);
MMS(:,2) = Index_MS_sorted(201:220,1);
HMS(:,1) = MS_sorted(276:295,1);
HMS(:,2) = Index_MS_sorted(276:295,1);
figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name',  'Event PLOT', ...
      'Position', [0.1 0.3 0.5 0.5]);
 AX1 = subplot(3,1,1);
 for ROI = 1 : 96
     Time_event = find(Event{LMS(1,2)}(:,ROI) > 0);
     Index_ROI = ROI.*ones(1,length(Time_event)); 
     plot(Time_event,Index_ROI,'b.');
     hold on
 end
xlabel('time(vol.)', ...
            'FontName', 'Arial');
ylabel('ROI', ...
            'FontName', 'Arial');
xlim([0,1200]);
ylim([0,96]);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize',18);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('MS=0.2394', ...
      'FontName','Arial', ...
      'FontSize', 20, ...
      'FontWeight', 'bold');

AX1 = subplot(3,1,2);
 for ROI = 1 : 96
     Time_event = find(Event{MMS(1,2)}(:,ROI) > 0);
     Index_ROI = ROI.*ones(1,length(Time_event)); 
     plot(Time_event,Index_ROI,'b.');
     hold on
 end
xlabel('time(vol.)', ...
            'FontName', 'Arial');
ylabel('ROI', ...
            'FontName', 'Arial');
xlim([0,1200]);
ylim([0,96]);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize',18);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('MS=0.4983', ...
      'FontName','Arial', ...
      'FontSize', 20, ...
      'FontWeight', 'bold');
 AX1 = subplot(3,1,3);
 for ROI = 1 : 96
     Time_event = find(Event{HMS(20,2)}(:,ROI) > 0);
     Index_ROI = ROI.*ones(1,length(Time_event)); 
     plot(Time_event,Index_ROI,'b.');
     hold on
 end
xlabel('time(vol.)', ...
            'FontName', 'Arial');
ylabel('ROI', ...
            'FontName', 'Arial');
xlim([0,1200]);
ylim([0,96]);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize',18);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('MS=0.6866', ...
      'FontName','Arial', ...
      'FontSize', 20, ...
      'FontWeight', 'bold');
%% avalanches analysis for LMS group
clc;clear;close all;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\HY_96_2');
Threshold = 1.4;
load('D:\supplyment_dataANDcode\c_MS_SE_static\MS_SE_HY96_RS.mat', 'syn')
[MS_sorted, Index_MS_sorted] = sort(syn); % rank the mean synchronize
LMS(:,1) = MS_sorted(1:20,1);
LMS(:,2) = Index_MS_sorted(1:20,1);
AvaS_LMS = [];
AvaD_LMS = []; 
BinE_LMS = [];
I_T = 0;
for I = 1:20
    S = LMS(I,2);
    I_T = I_T + 1;
    disp(['subject=',num2str(S)]);
    %% load the ROIs' signals and normalize the signals
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = rest_HY96_ROI;
    % z-normailiez 
    Signal_normalized = z_normalize(Signals_original);
    %% Avalanche analysis
    % avalanche define and ststistic
    BinSize=1;
    [Event{S,1}, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
    % the group-aggregate size, duration, events
    AvaS_LMS = [AvaS_LMS; AvaS{S,I_T}];
    AvaD_LMS = [AvaD_LMS; DurT{S,I_T}];
    BinE_LMS = [BinE_LMS; BinE{S,I_T}];
end
save(['Avalanches_LMS_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_LMS', 'AvaD_LMS', 'BinE_LMS', 'Event');
%% avalanches analysis for MMS group
clc;clear;close all;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\HY_96_2');
Threshold = 1.4;
load('D:\supplyment_dataANDcode\c_MS_SE_static\MS_SE_HY96_RS.mat', 'syn')
[MS_sorted, Index_MS_sorted] = sort(syn); % rank the mean synchroniz
MMS(:,1) = MS_sorted(201:220,1);
MMS(:,2) = Index_MS_sorted(201:220,1);
AvaS_MMS = [];
AvaD_MMS = []; 
BinE_MMS = [];
I_T = 0;
for I = 1:20
    S = MMS(I,2);
    I_T = I_T + 1;
    disp(['subject=',num2str(S)]);
    %% load the ROIs' signals and normalize the signals
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = rest_HY96_ROI;
    % z-normailiez 
    Signal_normalized = z_normalize(Signals_original);
    %% Avalanche analysis
    % avalanche define and ststistic
    BinSize=1;
    [Event{S,1}, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
    % the group-aggregate size, duration, events
    AvaS_MMS = [AvaS_MMS; AvaS{S,I_T}];
    AvaD_MMS = [AvaD_MMS; DurT{S,I_T}];
    BinE_MMS = [BinE_MMS; BinE{S,I_T}];
end
save(['Avalanches_MMS_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_MMS', 'AvaD_MMS', 'BinE_MMS', 'Event');
%% avalanches analysis for HMS group
clc;clear;close all;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\HY_96_2')
Threshold = 1.4;
load('D:\supplyment_dataANDcode\c_MS_SE_static\MS_SE_HY96_RS.mat', 'syn')
[MS_sorted, Index_MS_sorted] = sort(syn); % rank the mean synchroniz
HMS(:,1) = MS_sorted(276:295,1);
HMS(:,2) = Index_MS_sorted(276:295,1);
AvaS_HMS = [];
AvaD_HMS = []; 
BinE_HMS = [];
I_T = 0;
for I = 1:20
    S = HMS(I,2);
    I_T = I_T + 1;
    disp(['subject=',num2str(S)]);
    %% load the ROIs' signals and normalize the signals
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = rest_HY96_ROI;
    % z-normailiez 
    Signal_normalized = z_normalize(Signals_original);
    %% Avalanche analysis
    % avalanche define and ststistic
    BinSize=1;
    [Event{S,1}, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
    % the group-aggregate size, duration, events
    AvaS_HMS = [AvaS_HMS; AvaS{S,I_T}];
    AvaD_HMS = [AvaD_HMS; DurT{S,I_T}];
    BinE_HMS = [BinE_HMS; BinE{S,I_T}];
end
save(['Avalanches_HMS_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_HMS', 'AvaD_HMS', 'BinE_HMS', 'Event');
%%
clc;clear;close all
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\HY_96_2');
load('Avalanches_LMS_1.4.mat')
load('Avalanches_MMS_1.4.mat')
load('Avalanches_HMS_1.4.mat')
figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name',  'Power law analysis', ...
      'Position', [0.1 0.3 0.8 0.5]);
AX1 = subplot(1,3,1);
SizeRange_LMS = min(AvaS_LMS) : 1 : max(AvaS_LMS);
AvaSizeHist_LMS = hist(AvaS_LMS, SizeRange_LMS);
AvaSizeHist_LMS = AvaSizeHist_LMS/AvaSizeHist_LMS(1);
loglog(SizeRange_LMS, AvaSizeHist_LMS, ...
    'Color', [0, 0, 1], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
SizeRange_MMS = min(AvaS_MMS) : 1 : max(AvaS_MMS);
AvaSizeHist_MMS = hist(AvaS_MMS, SizeRange_MMS);
AvaSizeHist_MMS = AvaSizeHist_MMS/AvaSizeHist_MMS(1);
loglog(SizeRange_MMS, AvaSizeHist_MMS, ...
    'Color', [0, 1, 0], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
SizeRange_HMS = min(AvaS_HMS) : 1 : max(AvaS_HMS);
AvaSizeHist_HMS = hist(AvaS_HMS, SizeRange_HMS);
AvaSizeHist_HMS = AvaSizeHist_HMS/AvaSizeHist_HMS(1);
loglog(SizeRange_HMS, AvaSizeHist_HMS, ...
    'Color', [1, 0, 0], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
set(gca, ...
         'FontName', 'Arial', ...
         'FontSize',18, ...
         'XTick', [1,10,100,1000], ...
         'YTick', [0.00001,0.0001,0.001,0.01,0.1,1]);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off

AX1 = subplot(1,3,2);
DurRange_LMS = min(AvaD_LMS) : 1 : max(AvaD_LMS);
AvaDurHist_LMS = hist(AvaD_LMS, DurRange_LMS);
AvaDurHist_LMS = AvaDurHist_LMS/AvaDurHist_LMS(1);
loglog(DurRange_LMS, AvaDurHist_LMS, ...
    'Color', [0, 0, 1], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
DurRange_MMS = min(AvaD_MMS) : 1 : max(AvaD_MMS);
AvaDurHist_MMS = hist(AvaD_MMS, DurRange_MMS);
AvaDurHist_MMS = AvaDurHist_MMS/AvaDurHist_MMS(1);
loglog(DurRange_MMS, AvaDurHist_MMS, ...
    'Color', [0, 1, 0], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
DurRange_HMS = min(AvaD_HMS) : 1 : max(AvaD_HMS);
AvaDurHist_HMS = hist(AvaD_HMS, DurRange_HMS);
AvaDurHist_HMS = AvaDurHist_HMS/AvaDurHist_HMS(1);
loglog(DurRange_HMS, AvaDurHist_HMS, ...
    'Color', [1, 0, 0], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
hold on;
set(gca, ...
         'FontName', 'Arial', ...
         'FontSize',18, ...
         'XTick', [1,10,100,1000], ...
         'YTick', [0.00001,0.0001,0.001,0.01,0.1,1]);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off

AX1 = subplot(1,3,3);
[Lifetime_LMS, AverageSize_LMS]=lifetimeAverageSize(AvaD_LMS, AvaS_LMS);
loglog(Lifetime_LMS, AverageSize_LMS, ...
             'Color', [0,0,1], ...
             'LineStyle', 'none', ...
             'LineWidth', 2, ...
             'Marker', 'o', ...
             'MarkerSize',6);
hold on
[Lifetime_MMS, AverageSize_MMS]=lifetimeAverageSize(AvaD_MMS, AvaS_MMS);
loglog(Lifetime_MMS, AverageSize_MMS, ...
             'Color', [0,1,0], ...
             'LineStyle', 'none', ...
             'LineWidth', 2, ...
             'Marker', 'o', ...
             'MarkerSize',6);
hold on
[Lifetime_HMS, AverageSize_HMS]=lifetimeAverageSize(AvaD_HMS, AvaS_HMS);
loglog(Lifetime_HMS, AverageSize_HMS, ...
             'Color', [1,0,0], ...
             'LineStyle', 'none', ...
             'LineWidth', 2, ...
             'Marker', 'o', ...
             'MarkerSize',6);
set(gca, ...
         'FontName', 'Arial', ...
         'FontSize',18);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
%%
clc; clear; close all
load('Avalanches_LMS_1.4.mat')
addpath('D:\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25')
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_LMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_LMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end

[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_LMS, AvaS_LMS);
for Dura_min = 1 : 5
    for Dura_max = 7 : 20
        X = Dura_min
        Y = Dura_max-6
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_LMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_LMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('powerlaw_diff_fitarea_1.4_HY96LMS.mat');
%%
clc;clear;close all
load('powerlaw_diff_fitarea_1.4_HY96LMS.mat')
Alpha_r = reshape(Alpha,I*J,1);
P_alpha_r = reshape(P_alpha,I*J,1);
Smin2fit_r = reshape(Smin2fit,I*J,1);
Smax2fit_r = reshape(Smax2fit,I*J,1);
Alpha_rc = Alpha_r(find(P_alpha_r>0.2));
[Alpha_rc,Alpha_rc_index] = sort(Alpha_rc);
P_alpha_rc = P_alpha_r(find(P_alpha_r>0.2));
P_alpha_rc = P_alpha_rc(Alpha_rc_index);
Smin2fit_rc = Smin2fit_r(find(P_alpha_r>0.2));
Smin2fit_rc = Smin2fit_rc(Alpha_rc_index);
Smax2fit_rc =Smax2fit_r(find(P_alpha_r>0.2));
Smax2fit_rc = Smax2fit_rc(Alpha_rc_index);

Tau_r = reshape(Tau,X*Y,1);
P_tau_r = reshape(P_tau,X*Y,1);
Tmin2fit_r = reshape(Tmin2fit,X*Y,1);
Tmax2fit_r = reshape(Tmax2fit,X*Y,1);
Gamma_fit_r = reshape(Gamma_fit,X*Y,1);
Gamma_fit_a_r = reshape(Gamma_fit_a,X*Y,1);
Tau_rc = Tau_r(find(P_tau_r>0.2));
[Tau_rc,Tau_rc_index] = sort(Tau_rc);
P_tau_rc = P_tau_r(find(P_tau_r>0.2));
P_tau_rc = P_tau_rc(Tau_rc_index);
Tmin2fit_rc = Tmin2fit_r(find(P_tau_r>0.2));
Tmin2fit_rc = Tmin2fit_rc(Tau_rc_index);
Tmax2fit_rc = Tmax2fit_r(find(P_tau_r>0.2));
Tmax2fit_rc = Tmax2fit_rc(Tau_rc_index);
Gamma_fit_rc = Gamma_fit_r(find(P_tau_r>0.2));
Gamma_fit_rc = Gamma_fit_rc(Tau_rc_index);
Gamma_fit_a_rc = Gamma_fit_a_r(find(P_tau_r>0.2));
Gamma_fit_a_rc = Gamma_fit_a_rc(Tau_rc_index);

for M = 1 : length(Alpha_rc)
    for N = 1 : length(Tau_rc) 
        Gamma_alpha_tau(M,N) = (Tau_rc(N)-1) / (Alpha_rc(M)-1);
        Gamma_fit_1(M,N) = Gamma_fit_rc(N);
    end
end
DLAT = abs(Gamma_alpha_tau-Gamma_fit_1);

AX1 = subplot(3,2,1);
imagesc([30,60],[1,10],Alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\alpha$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,2);
imagesc([30,60],[1,10],P_alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,4);
imagesc([7 20],[1 5],Tau);
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\tau$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,5);
imagesc([7 20],[1 5],P_tau)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,6);
imagesc([7 20],[1 5],Gamma_fit)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\gamma$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');

AX1 = subplot(3,2,5);
imagesc(Gamma_alpha_tau)
xticks([1 5 10 15 20 25])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(5)),num2str(Tau_rc(10)),...
    num2str(Tau_rc(15)),num2str(Tau_rc(20)),num2str(Tau_rc(25))});
yticks([1 25 50 75 100 125 150 175])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(25)),num2str(Alpha_rc(50)),...
    num2str(Alpha_rc(75)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(125)),num2str(Alpha_rc(150)),num2str(Alpha_rc(175))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\frac{1-\tau}{1-\alpha}$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,6);
N_delta_01 = 0;
for X = 1 : 178
    for Y = 1 : 25
        if DLAT(X,Y) <= 0.1
            N_delta_01 = N_delta_01 + 1;
           Delta(X,Y) = DLAT(X,Y);
        end
        if DLAT(X,Y) > 0.1
            Delta(X,Y) = 0.15;
        end
    end
end
imagesc(Delta)
xticks([1 5 10 15 20 25])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(5)),num2str(Tau_rc(10)),...
    num2str(Tau_rc(15)),num2str(Tau_rc(20)),num2str(Tau_rc(25))});
yticks([1 25 50 75 100 125 150 175])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(25)),num2str(Alpha_rc(50)),...
    num2str(Alpha_rc(75)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(125)),num2str(Alpha_rc(150)),num2str(Alpha_rc(175))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$|\frac{1-\tau}{1-\alpha}-\gamma|$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
%%
clc; clear; close all
load('Avalanches_MMS_1.4.mat')
SizeFitRange = [3,30];
DurationFitRange = [2,9];
[Avalanche, Criterion, BP_GROUP] = CriticalAnalysis_powerlaw(AvaD_MMS, AvaS_MMS, BinE_MMS, SizeFitRange, DurationFitRange);
addpath('D:\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25')
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_MMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_MMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end

[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_MMS, AvaS_MMS);
for Dura_min = 1 : 5
    for Dura_max = 7 : 20
        X = Dura_min
        Y = Dura_max-6
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_MMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_MMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('powerlaw_diff_fitarea_1.4_HY96MMS.mat');
%%
clc;clear;close all
load('powerlaw_diff_fitarea_1.4_HY96MMS.mat')
P_thr = 0.2;
Alpha_r = reshape(Alpha,I*J,1);
P_alpha_r = reshape(P_alpha,I*J,1);
Smin2fit_r = reshape(Smin2fit,I*J,1);
Smax2fit_r = reshape(Smax2fit,I*J,1);
Alpha_rc = Alpha_r(find(P_alpha_r>P_thr));
[Alpha_rc,Alpha_rc_index] = sort(Alpha_rc);
P_alpha_rc = P_alpha_r(find(P_alpha_r>P_thr));
P_alpha_rc = P_alpha_rc(Alpha_rc_index);
Smin2fit_rc = Smin2fit_r(find(P_alpha_r>P_thr));
Smin2fit_rc = Smin2fit_rc(Alpha_rc_index);
Smax2fit_rc =Smax2fit_r(find(P_alpha_r>P_thr));
Smax2fit_rc = Smax2fit_rc(Alpha_rc_index);

Tau_r = reshape(Tau,X*Y,1);
P_tau_r = reshape(P_tau,X*Y,1);
Tmin2fit_r = reshape(Tmin2fit,X*Y,1);
Tmax2fit_r = reshape(Tmax2fit,X*Y,1);
Gamma_fit_r = reshape(Gamma_fit,X*Y,1);
Gamma_fit_a_r = reshape(Gamma_fit_a,X*Y,1);
Tau_rc = Tau_r(find(P_tau_r>0.2));
[Tau_rc,Tau_rc_index] = sort(Tau_rc);
P_tau_rc = P_tau_r(find(P_tau_r>P_thr));
P_tau_rc = P_tau_rc(Tau_rc_index);
Tmin2fit_rc = Tmin2fit_r(find(P_tau_r>P_thr));
Tmin2fit_rc = Tmin2fit_rc(Tau_rc_index);
Tmax2fit_rc = Tmax2fit_r(find(P_tau_r>P_thr));
Tmax2fit_rc = Tmax2fit_rc(Tau_rc_index);
Gamma_fit_rc = Gamma_fit_r(find(P_tau_r>P_thr));
Gamma_fit_rc = Gamma_fit_rc(Tau_rc_index);
Gamma_fit_a_rc = Gamma_fit_a_r(find(P_tau_r>P_thr));
Gamma_fit_a_rc = Gamma_fit_a_rc(Tau_rc_index);

for M = 1 : length(Alpha_rc)
    for N = 1 : length(Tau_rc) 
        Gamma_alpha_tau(M,N) = (Tau_rc(N)-1) / (Alpha_rc(M)-1);
        Gamma_fit_1(M,N) = Gamma_fit_rc(N);
    end
end
DLAT = abs(Gamma_alpha_tau-Gamma_fit_1);

AX1 = subplot(3,2,1);
imagesc([30,60],[1,10],Alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\alpha$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,2);
imagesc([30,60],[1,10],P_alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,4);
imagesc([7 20],[1 5],Tau);
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\tau$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,5);
imagesc([7 20],[1 5],P_tau)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,6);
imagesc([7 20],[1 5],Gamma_fit)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\gamma$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');

AX1 = subplot(3,2,5);
imagesc(Gamma_alpha_tau)
xticks([1 5 15 25 35 45])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(5)),num2str(Tau_rc(15)),...
    num2str(Tau_rc(25)),num2str(Tau_rc(35)),num2str(Tau_rc(45))})
yticks([1 50 100 150 200 250 300])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(150)),num2str(Alpha_rc(200)),num2str(Alpha_rc(250)),num2str(Alpha_rc(300))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\frac{1-\tau}{1-\alpha}$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,6);
N_delta_01 = 0;
for X = 1 : 310
    for Y = 1 : 47
        if DLAT(X,Y) <= 0.1
           Delta(X,Y) = DLAT(X,Y);
           N_delta_01 = N_delta_01 +1;
        end
        if DLAT(X,Y) > 0.1
            Delta(X,Y) = 0.15;
        end
    end
end
imagesc(Delta)
xticks([1 5 15 25 35 45])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(5)),num2str(Tau_rc(15)),...
    num2str(Tau_rc(25)),num2str(Tau_rc(35)),num2str(Tau_rc(45))})
yticks([1 50 100 150 200 250 300])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(150)),num2str(Alpha_rc(200)),num2str(Alpha_rc(250)),num2str(Alpha_rc(300))});
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$|\frac{1-\tau}{1-\alpha}-\gamma|$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
%%
clc; clear; close all
load('Avalanches_HMS_1.4.mat')
SizeFitRange = [3,30];
DurationFitRange = [2,9];
addpath('D:\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25')
for Size_min = 1 : 10
    for Size_max = 30 : 60
        I = Size_min
        J = Size_max-29;
        [Alpha(I, J), Smin2fit(I, J), Smax2fit(I, J)] = plmle(AvaS_HMS, 'xmin', Size_min, 'xmax', Size_max);
        [P_alpha(I, J)] = pvcalc(AvaS_HMS, Alpha(I, J), 'xmin',Smin2fit(I, J),'xmax',Smax2fit(I, J));
    end
end

[Lifetime, AverageSize]=lifetimeAverageSize(AvaD_HMS, AvaS_HMS);
for Dura_min = 1 : 5
    for Dura_max = 7 : 20
        X = Dura_min
        Y = Dura_max-6;
        [Tau(X, Y), Tmin2fit(X, Y), Tmax2fit(X, Y)] = ...
            plmle(AvaD_HMS,'xmin', Dura_min, 'xmax', Dura_max);
        
        [P_tau(X, Y)] = ...
            pvcalc(AvaD_HMS, Tau(X, Y), 'xmin',Tmin2fit(X, Y),'xmax',Tmax2fit(X, Y));
        
        [Gamma_nlsfit] = ...
            createFit_pl_nonlinearLS(Lifetime(Tmin2fit(X, Y):Tmax2fit(X, Y)), AverageSize(Tmin2fit(X, Y):Tmax2fit(X, Y)));
        
        Gamma_fit(X, Y) = Gamma_nlsfit.b;
        Gamma_fit_a(X, Y) = Gamma_nlsfit.a;
    end
end
save('powerlaw_diff_fitarea_1.4_HY96HMS.mat')
%%
clc;clear;close all
load('powerlaw_diff_fitarea_1.4_HY96HMS.mat')
Alpha_r = reshape(Alpha,I*J,1);
P_alpha_r = reshape(P_alpha,I*J,1);
Smin2fit_r = reshape(Smin2fit,I*J,1);
Smax2fit_r = reshape(Smax2fit,I*J,1);
Alpha_rc = Alpha_r(find(P_alpha_r>0.2));
[Alpha_rc,Alpha_rc_index] = sort(Alpha_rc);
P_alpha_rc = P_alpha_r(find(P_alpha_r>0.2));
P_alpha_rc = P_alpha_rc(Alpha_rc_index);
Smin2fit_rc = Smin2fit_r(find(P_alpha_r>0.2));
Smin2fit_rc = Smin2fit_rc(Alpha_rc_index);
Smax2fit_rc =Smax2fit_r(find(P_alpha_r>0.2));
Smax2fit_rc = Smax2fit_rc(Alpha_rc_index);

Tau_r = reshape(Tau,X*Y,1);
P_tau_r = reshape(P_tau,X*Y,1);
Tmin2fit_r = reshape(Tmin2fit,X*Y,1);
Tmax2fit_r = reshape(Tmax2fit,X*Y,1);
Gamma_fit_r = reshape(Gamma_fit,X*Y,1);
Gamma_fit_a_r = reshape(Gamma_fit_a,X*Y,1);
Tau_rc = Tau_r(find(P_tau_r>0.2));
[Tau_rc,Tau_rc_index] = sort(Tau_rc);
P_tau_rc = P_tau_r(find(P_tau_r>0.2));
P_tau_rc = P_tau_rc(Tau_rc_index);
Tmin2fit_rc = Tmin2fit_r(find(P_tau_r>0.2));
Tmin2fit_rc = Tmin2fit_rc(Tau_rc_index);
Tmax2fit_rc = Tmax2fit_r(find(P_tau_r>0.2));
Tmax2fit_rc = Tmax2fit_rc(Tau_rc_index);
Gamma_fit_rc = Gamma_fit_r(find(P_tau_r>0.2));
Gamma_fit_rc = Gamma_fit_rc(Tau_rc_index);
Gamma_fit_a_rc = Gamma_fit_a_r(find(P_tau_r>0.2));
Gamma_fit_a_rc = Gamma_fit_a_rc(Tau_rc_index);

for M = 1 : length(Alpha_rc)
    for N = 1 : length(Tau_rc) 
        Gamma_alpha_tau(M,N) = (Tau_rc(N)-1) / (Alpha_rc(M)-1);
        Gamma_fit_1(M,N) = Gamma_fit_rc(N);
    end
end
DLAT = abs(Gamma_alpha_tau-Gamma_fit_1);

AX1 = subplot(3,2,1);
imagesc([30,60],[1,10],Alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\alpha$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,2);
imagesc([30,60],[1,10],P_alpha)
xlabel('maximal size S', ...
       'FontName', 'Arial');
ylabel('minimal size S', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,4);
imagesc([7 20],[1 5],Tau);
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\tau$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,5);
imagesc([7 20],[1 5],P_tau)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('P', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,3,6);
imagesc([7 20],[1 5],Gamma_fit)
xlabel('maximal duration T', ...
       'FontName', 'Arial');
ylabel('minimal duration T', ...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\gamma$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');

AX1 = subplot(3,2,5);
Gamma_alpha_tau_1 = Gamma_alpha_tau;
for X = 1 : 310
    for Y = 1 :61
        if Gamma_alpha_tau(X,Y) >= max(max(Gamma_fit+5))
           Gamma_alpha_tau_1(X,Y) = max(max(Gamma_fit+5));
        end
    end
end
imagesc(Gamma_alpha_tau_1)
xticks([1 10 20 30 40 50 60])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(10)),num2str(Tau_rc(20)),...
    num2str(Tau_rc(30)),num2str(Tau_rc(40)),num2str(Tau_rc(50)),num2str(Tau_rc(60))})
yticks([1 50 100 150 200 250 300])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(150)),num2str(Alpha_rc(200)),num2str(Alpha_rc(250)),num2str(Alpha_rc(300))})
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$\frac{1-\tau}{1-\alpha}$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');
  
AX1 = subplot(3,2,6);
N_delta_01 = 0;
for X = 1 : 310
    for Y = 1 :61
        if DLAT(X,Y) <= 0.1
           N_delta_01 = N_delta_01 + 1;
           Delta(X,Y) = DLAT(X,Y);
        end
        if DLAT(X,Y) > 0.1
            Delta(X,Y) = 0.15;
        end
    end
end
imagesc(Delta)
xticks([1 10 20 30 40 50 60])
xticklabels({num2str(Tau_rc(1)),num2str(Tau_rc(10)),num2str(Tau_rc(20)),...
    num2str(Tau_rc(30)),num2str(Tau_rc(40)),num2str(Tau_rc(50)),num2str(Tau_rc(60))})
yticks([1 50 100 150 200 250 300])
yticklabels({num2str(Alpha_rc(1)),num2str(Alpha_rc(50)),num2str(Alpha_rc(100)),...
    num2str(Alpha_rc(150)),num2str(Alpha_rc(200)),num2str(Alpha_rc(250)),num2str(Alpha_rc(300))})
xlabel('$\tau$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
ylabel('$\alpha$', ...
      'Interpreter','latex',...
       'FontName', 'Arial');
set(gca, ...
    'FontName', 'Arial', ...
     'FontSize',13);
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off
title('$|\frac{1-\tau}{1-\alpha}-\gamma|$', ...
      'Interpreter','latex',...
      'FontName','Arial', ...
      'FontSize', 18, ...
      'FontWeight', 'bold');