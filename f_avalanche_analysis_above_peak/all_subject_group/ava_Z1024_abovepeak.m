%% The first part: avalanche define
clc;clear;close all
addpath('D:\the_code_and_results_sort\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25');
I_T = 0;
for Threshold = 1 : 0.1 : 2.9
    I_T = I_T + 1;
    AvaS_AllSubject = [];
    AvaD_AllSubject = []; 
    BinE_AllSubject = [];
    for S=1:295 % subject
        disp(['subject=',num2str(S)]);
       %% load the ROIs' signals and normalize the signals
        load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(S),'.mat']); %load the ROI signals
        Signals_original = ROI_1024_RS;
        % z-normailiez 
        Signal_normalized = z_normalize(Signals_original);
       %% Avalanche analysis
        % avalanche define and ststistic
        BinSize=1;
        [Event, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
        % the group-aggregate size, duration, events
        AvaS_AllSubject = [AvaS_AllSubject; AvaS{S,I_T}];
        AvaD_AllSubject = [AvaD_AllSubject; DurT{S,I_T}];
        BinE_AllSubject = [BinE_AllSubject; BinE{S,I_T}];
        jj=0;
        for ROI = 1 : 1024
            for T = 1 : 1200
                if Event(T,ROI) == 1
                    jj = jj + 1;
                    Event_serise{S,1}(jj,1) = ROI;
                    Event_serise{S,1}(jj,2) = T;
                end
            end
        end
    end
    BP_subject(I_T) = branchingParameter(BinE_AllSubject);
    save(['Z_1024_1\Avalanches_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_AllSubject', 'AvaD_AllSubject', 'BinE_AllSubject');
    save(['Z_1024_1\Z1024_event_', num2str(Threshold), '.mat'],'Event_serise');
end
save('Z_1024_1\Avalanches_define.mat', 'BinE', 'AvaS', 'DurT', 'AvaSH', 'AvaDH', 'BP', 'BP_other');
save('Z_1024_1\BP.mat','BP_subject');
%% The first part: avalanche define
clc;clear;close all
addpath('D:\the_code_and_results_sort\supplyment_dataANDcode\z_function\CritAnalysisSoftwarePackage2016-04-25');
I_T = 0;
for Threshold = 2 : 0.1 : 2
    I_T = I_T + 1;
    AvaS_AllSubject = [];
    AvaD_AllSubject = []; 
    BinE_AllSubject = [];
    for S=1:295 % subject
        disp(['subject=',num2str(S)]);
       %% load the ROIs' signals and normalize the signals
        load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(S),'.mat']); %load the ROI signals
        Signals_original = ROI_1024_RS;
        % z-normailiez 
        Signal_normalized = z_normalize(Signals_original);
       %% Avalanche analysis
        % avalanche define and ststistic
        BinSize=1;
        [Event, BinE{S,I_T}, AvaS{S,I_T}, DurT{S,I_T}, AvaSH{S,I_T},AvaDH{S,I_T}, BP(S,I_T), BP_other(S,I_T,:)] =...
                   avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
        % the group-aggregate size, duration, events
        AvaS_AllSubject = [AvaS_AllSubject; AvaS{S,I_T}];
        AvaD_AllSubject = [AvaD_AllSubject; DurT{S,I_T}];
        BinE_AllSubject = [BinE_AllSubject; BinE{S,I_T}];
        jj=0;
        for ROI = 1 : 1024
            for T = 1 : 1200
                if Event(T,ROI) == 1
                    jj = jj + 1;
                    Event_serise{S,1}(jj,1) = ROI;
                    Event_serise{S,1}(jj,2) = T;
                end
            end
        end
    end
    BP_subject(I_T) = branchingParameter(BinE_AllSubject);
    save(['Z_1024\Avalanches_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_AllSubject', 'AvaD_AllSubject', 'BinE_AllSubject');
    save(['Z_1024\Z1024_event_', num2str(Threshold), '.mat'],'Event_serise');
end
save('Z_1024\Avalanches_define.mat', 'BinE', 'AvaS', 'DurT', 'AvaSH', 'AvaDH', 'BP', 'BP_other');
save('Z_1024\BP.mat','BP_subject');
%% avalanche analysis for LMS group in z1024 
clc;clear;close all;
Threshold = 2.6;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\Z1024')
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
    load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = ROI_1024_RS;
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
%% avalanche analysis for MMS group in z1024 
clc;clear;close all;
Threshold = 2.6;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\Z1024')
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
     load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = ROI_1024_RS;
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
%% avalanche analysis for HMS group in z1024 
clc;clear;close all;
Threshold = 2.6;
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\Z1024')
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
     load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\Z_1024\sub',num2str(S),'.mat']); %load the ROI signals
    Signals_original = ROI_1024_RS;
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
%% power law analysis for LMS, MMS and HMS
clc;clear;close all
cd('D:\supplyment_dataANDcode\b_avalanche_analysis_grouply\Z1024');
load('Avalanches_LMS_2.6.mat')
load('Avalanches_MMS_2.6.mat')
load('Avalanches_HMS_2.6.mat')
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
