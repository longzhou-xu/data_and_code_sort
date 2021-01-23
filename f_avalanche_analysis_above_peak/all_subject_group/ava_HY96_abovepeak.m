%% The first part: avalanche define
clc;clear;close all
addpath(['D:\the_code_and_results_sort\supplyment_dataANDcode\',...
    'z_function\CritAnalysisSoftwarePackage2016-04-25']);
I_T = 0;
THD = 1:0.1:2;
for Threshold = 1 : 0.1 : 2 
    I_T = I_T + 1;
    AvaS_AllSubject = [];
    AvaD_AllSubject = []; 
    BinE_AllSubject = [];
    for S=1:295 % subject
        disp(['subject=',num2str(S)]);
       %% load the ROIs' signals and normalize the signals
        load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
            'a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
            num2str(S),'.mat']); %load the ROI signals
        Signals_original = rest_HY96_ROI;
        % z-normailiez 
        Signal_normalized = z_normalize(Signals_original);
       %% Avalanche analysis
        % avalanche define and ststistic
        BinSize=1;
        [Event, BinE, AvaS, DurT, AvaSH,...
            AvaDH, BP(S,1), BP_other(S,:)] =...
            avalancheAnalysisBOLD(Signal_normalized,Threshold, BinSize);
        % the group-aggregate size, duration, events
        AvaS_AllSubject = [AvaS_AllSubject; AvaS];
        AvaD_AllSubject = [AvaD_AllSubject; DurT];
        BinE_AllSubject = [BinE_AllSubject; BinE];
        %event serise
        clear AvaS DurT BinE
        jj=0;
        for ROI = 1 : 96
            for T = 1 : 1200
                if Event(T,ROI) == 1
                    jj = jj + 1;
                    Event_serise{S,1}(jj,1) = ROI;
                    Event_serise{S,1}(jj,2) = T;
                end
            end
        end
        % event cv
        Event_sum(1:1200,S) = sum(Event,2);
        Event_rate(1:1200,S) = Event_sum(1:1200,S)/96;
        Event_CV(S) = std(Event_rate(1:1200,S))/mean(Event_rate(1:1200,S));
        clear Event 
    end
    save(['HY_96\Ava_diff_thd\Avalanches_', num2str(Threshold), '.mat'], ...
        'Threshold', 'AvaS_AllSubject', 'AvaD_AllSubject', 'BinE_AllSubject',...
        'Event_serise','BP','BP_other','Event_sum','Event_rate','Event_CV');
    BP_GROUP(I_T) = branchingParameter(BinE_AllSubject);
end
save('HY_96\BP_group.mat','BP_GROUP','THD');
%% power law analysis
clc; clear; close all
load('HY_96\Ava_diff_thd\Avalanches_2.mat')
SizeFitRange = [3,30];
DurationFitRange = [3,9];
[Avalanche, Criterion, BP_GROUP] = CriticalAnalysis_powerlaw(...
    AvaD_AllSubject, AvaS_AllSubject, BinE_AllSubject, SizeFitRange,...
    DurationFitRange);
% plot the bp
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_HY96_RS.mat'], 'syn')
figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name',  'Power law analysis', ...
      'Position', [0.1 0.3 0.5 0.5]);
AX1 = subplot(1,1,1);
plot(syn, BP(:,1) ,...
    'Color', [0.50, 0.50, 0.50], ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize',6);
[R,P] = corr(syn, BP(:,1));
xlabel('MS');
ylabel('BP')
grid off; 
AX1.LineWidth = 2; 
box(AX1,'off'); 
hold off

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
