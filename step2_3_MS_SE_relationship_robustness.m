%% the robust relationship between the MS and SE
% setting
clc;clear;close all;
figure('Color', 'w', 'Units', 'Normalized', ...
    'Name', ['The robustness of inverted-U curve',...
    ' between MS and SE against choice of parcellation'], ...
    'Position', [0.1 0.1 0.8 0.8]);

%% select the LMS-MMS-HMS group
load('step_1\subject.mat')
load('step_2\MS_SE_HY96_RS.mat', 'MS');
% sort the MS in ascending order
[MS_s, MS_s_i] = sort(MS);
% LMS group
LMS(:,1) = subject(MS_s_i(1:20));
LMS(:,2) = MS_s(1:20);
LMS(:,3) = MS_s_i(1:20);
% MMS group
MMS(:,1) = subject(MS_s_i(200:219));
MMS(:,2) = MS_s(200:219);
MMS(:,3) = MS_s_i(200:219);
% HMS group
HMS(:,1) = subject(MS_s_i(276:295));
HMS(:,2) = MS_s(276:295);
HMS(:,3) = MS_s_i(276:295);
save('step_2\HMS_MMS_LMS_subject.mat', 'LMS', 'MMS', 'HMS');

%% show the relationship between MS and SE
%setting
clc;clear;close all;

% HY-96 
load('step_2\MS_SE_HY96_RS.mat') % load data
SE_HY96 = SE; MS_HY96 = MS;  
clear MS SE KOP
load('step_2\HMS_MMS_LMS_subject.mat') % load the MMS, LMS and HMS subject index
AX1 = subplot(2,3,1);
plot(MS_HY96, SE_HY96, ...
    'Color',[0.5,0.5,0.5],...
    'Marker','o', 'MarkerSize',6,...
    'LineWidth',2, 'LineStyle','none'); % plot the MS and SE
hold on
f2 = plot(MS_HY96(LMS(:,3)), SE_HY96(LMS(:,3)),...
    'Color',[0,0,1],...
    'Marker','o', 'MarkerSize',6,...
    'LineWidth',2, 'LineStyle','none');
f3 = plot(MS_HY96(MMS(:,3)),SE_HY96(MMS(:,3)),...
    'Color',[0,1,0],...
    'Marker','o', 'MarkerSize',6,...
    'LineWidth',2, 'LineStyle','none');
% marker the HMS
f4 = plot(MS_HY96(HMS(:,3)), SE_HY96(HMS(:,3)),...
    'Color',[1,0,0],...
    'Marker','o', 'MarkerSize',6,...
    'LineWidth',2, 'LineStyle','none');% marker the HMS
MS_SE_fit = createFit_poly2(MS_HY96,SE_HY96);
p1 = MS_SE_fit.p1; 
p2 = MS_SE_fit.p2; 
p3 = MS_SE_fit.p3;
clear x y; 
x = min(MS_HY96)-0.01:0.01:max(MS_HY96)+0.02; 
y = p1.*(x).^2+p2.*x+p3;
f1 = plot(x, y, ...
    'Color', [1,0,0], ...
    'LineStyle', '--', 'LineWidth', 2.5);
xlim([min(x)-0.05,max(x)+0.05]);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
legend([f1,f2,f3,f4],{'$fit$','LMS','MMS','HMS'},...
    'Interpreter','latex',...
    'Location','northeast',...
    'FontName','Arial', 'FontSize',12,...
    'EdgeColor',[1,1,1])
legend('boxoff')
xlabel('Atlas_9_6 MS <r>',...
    'FontName','Arial',...
    'FontSize',15);
ylabel('Atlas_9_6 SE H(r)',...
    'FontName','Arial',...
    'FontSize',15);
box(AX1,'off'); 
AX1.LineWidth=2; 
grid off; 
hold off;
title('a', ...
    'FontName', 'Arial',...
    'FontSize',24, ...
    'FontWeight', 'bold',...
    'units','normalized',...
    'position',[-1/18,1+1/10],...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
%%
load('MS_SE_BA246_RS.mat')
SE_BA246=synE; 
MS_BA246=syn; 
clear syn synE subject KOP
AX1=subplot(2,3,2);
plot(MS_BA246,SE_BA246,...
    'Color',[0.5,0.5,0.5],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
hold on
f2 = plot(MS_BA246(subcri_subject_number(:,3)),...
    SE_BA246(subcri_subject_number(:,3)),...
    'Color',[0,0,1],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
f3 = plot(MS_BA246(cri_subject_number(:,3)),...
    SE_BA246(cri_subject_number(:,3)),...
    'Color',[0,1,0],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
f4 = plot(MS_BA246(suppercri_subject_number(:,3)),...
    SE_BA246(suppercri_subject_number(:,3)),...
    'Color',[1,0,0],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
MS_SE_fit=createFit_poly2(MS_BA246,SE_BA246);
p1=MS_SE_fit.p1; 
p2=MS_SE_fit.p2; 
p3=MS_SE_fit.p3;
clear x y; 
x=min(MS_BA246)-0.01:0.01:max(MS_BA246)+0.02; 
y=p1.*(x).^2+p2.*x+p3;
f1=plot(x,y,'Color',[1,0,0],'LineStyle','--','LineWidth',2.5);
xlim([min(x)-0.05,max(x)+0.05]);
set(gca,...
    'FontName','Arial',...
    'FontSize',12);
legend([f1,f2,f3,f4],{'$fit$','LMS','MMS','HMS'},...
    'Interpreter','latex',...
    'Location','southeast',...
    'FontName','Arial',...
    'FontSize',12,...
    'EdgeColor',[1,1,1])
legend('boxoff')
xlabel('Atlas_2_4_6 MS <r>',...
    'FontName','Arial',...
    'FontSize',15);
ylabel('Atlas_2_4_6 SE H(r)',...
    'FontName','Arial',...
    'FontSize',15);
box(AX1,'off'); 
AX1.LineWidth = 2; 
grid off; 
hold off;
title('b',...
    'FontName','Arial',...
    'FontSize',24, ...
    'FontWeight', 'bold',...
    'units','normalized',...
    'position',[-1/18,1+1/10],...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
%%
load('MS_SE_Z1024_RS.mat')
SE_Z1024 = synE; 
MS_Z1024 = syn;  
clear syn synE subject KOP
AX1 = subplot(2,3,3);
plot(MS_Z1024,SE_Z1024,...
    'Color',[0.5,0.5,0.5],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
hold on
f2 = plot(MS_Z1024(subcri_subject_number(:,3)),...
    SE_Z1024(subcri_subject_number(:,3)),...
    'Color',[0,0,1],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
f3 = plot(MS_Z1024(cri_subject_number(:,3)),....
    SE_Z1024(cri_subject_number(:,3)),...
    'Color',[0,1,0],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
f4 = plot(MS_Z1024(suppercri_subject_number(:,3)),...
    SE_Z1024(suppercri_subject_number(:,3)),...
    'Color',[1,0,0],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
MS_SE_fit=createFit_poly2(MS_Z1024,SE_Z1024);
p1 = MS_SE_fit.p1; 
p2 = MS_SE_fit.p2; 
p3 = MS_SE_fit.p3;
clear x y; 
x = min(MS_Z1024)-0.01:0.01:max(MS_Z1024)+0.02; 
y = p1.*(x).^2+p2.*x+p3;
f1 = plot(x,y,...
    'Color',[1,0,0],...
    'LineStyle','--',...
    'LineWidth',2.5);
xlim([min(x)-0.05,max(x)+0.05]);
set(gca,...
    'FontName','Arial',...
    'FontSize',12);
legend([f1,f2,f3,f4],{'$fit$','LMS','MMS','HMS'},...
    'Interpreter','latex',...
    'Location','southeast',...
    'FontName','Arial',...
    'FontSize',12,...
    'EdgeColor',[1,1,1])
legend('boxoff')
xlabel('Atlas_1_0_2_4 MS <r>',...
    'FontName','Arial',...
    'FontSize',15);
ylabel('Atlas_1_0_2_4 SE H(r)',...
    'FontName','Arial',...
    'FontSize',15);
box(AX1,'off'); 
AX1.LineWidth=2; 
grid off; 
hold off;
title('c',...
    'FontName','Arial',...
    'FontSize',24, ...
    'FontWeight', 'bold',...
    'units','normalized',...
    'position',[-1/18,1+1/10],...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
%%
AX1 = subplot(2,3,4);
plot(SE_HY96,SE_BA246,...
    'Color',[0.5,0.5,0.5],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
hold on
fit_SE96_SE246 = createFit_poly1(SE_HY96,SE_BA246);
[R_SE96_SE246,P_SE96_SE246] = corr(SE_HY96,SE_BA246);
p1 = fit_SE96_SE246.p1; 
p2 = fit_SE96_SE246.p2; 
fitX = min(SE_HY96)-0.1:0.1:max(SE_HY96)+0.1;
fitY = p1 .* fitX + p2;
f1 = plot(fitX, fitY, ...
    'Color', [1,0,0], ...
    'LineStyle', '--', ...
    'Linewidth', 2.5, ...
    'Marker', 'none', ...
    'MarkerSize', 2.5);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 12)
lgd = legend(f1, '$fit$', ...
    'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 12;
lgd.Location = 'northwest';
text(0.80, 0.35, ['$R=', num2str(round(R_SE96_SE246,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
text(0.80, 0.2, ['$p=', num2str(round(P_SE96_SE246,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
xlabel('Atlas_9_6 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
ylabel('Atlas_2_4_6 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
title('d', ...
    'FontName', 'Arial', ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'units','normalized',...
    'position', [-1/18,1+1/10], ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off
%%
AX1 = subplot(2,3,5);
SE_BA246(6)=[];
SE_HY96(6)=[];
SE_Z1024(6)=[];
MS_Z1024(6)=[];
plot(SE_HY96,SE_Z1024,...
    'Color',[0.5,0.5,0.5],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
hold on
fit_SE96_SE1024 = createFit_poly1(SE_HY96,SE_Z1024);
[R_SE96_SE1024,P_SE96_SE1024] = corr(SE_HY96,SE_Z1024);
p1 = fit_SE96_SE1024.p1; 
p2 = fit_SE96_SE1024.p2; 
fitX = min(SE_HY96)-0.1:0.1:max(SE_HY96)+0.1;
fitY = p1 .* fitX + p2;
f1 = plot(fitX, fitY, ...
    'Color', [1,0,0], ...
    'LineStyle', '--', ...
    'Linewidth', 2.5, ...
    'Marker', 'none', ...
    'MarkerSize', 2.5);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 12)
lgd = legend(f1, '$fit$', ...
    'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 12;
lgd.Location = 'northwest';
text(0.80, 0.35, ['$R=', num2str(round(R_SE96_SE1024,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
text(0.80, 0.2, ['$p=', num2str(round(P_SE96_SE1024,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
xlabel('Atlas_9_6 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
ylabel('Atlas_1_0_2_4 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
title('e', ...
    'FontName', 'Arial', ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'units','normalized',...
    'position', [-1/18,1+1/10], ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off
%%
AX1 = subplot(2,3,6);
plot(SE_BA246,SE_Z1024,...
    'Color',[0.5,0.5,0.5],...
    'Marker','o',...
    'MarkerSize',6,...
    'LineWidth',2,...
    'LineStyle','none');
hold on
fit_SE246_SE1024 = createFit_poly1(SE_BA246,SE_Z1024);
[R_SE246_SE1024,P_SE246_SE1024] = corr(SE_BA246,SE_Z1024);
p1 = fit_SE246_SE1024.p1; 
p2 = fit_SE246_SE1024.p2; 
fitX = min(SE_BA246)-0.1:0.1:max(SE_BA246)+0.1;
fitY = p1 .* fitX + p2;
f1 = plot(fitX, fitY, ...
    'Color', [1,0,0], ...
    'LineStyle', '--', ...
    'Linewidth', 2.5, ...
    'Marker', 'none', ...
    'MarkerSize', 2.5);
set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 12)
lgd = legend(f1, '$fit$', ...
    'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 12;
lgd.Location = 'northwest';
text(0.80, 0.35, ['$R=', num2str(round(R_SE246_SE1024,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
text(0.80, 0.2, ['$p=', num2str(round(P_SE246_SE1024,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
xlabel('Atlas_2_4_6 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
ylabel('Atlas_1_0_2_4 SE', ...
    'FontName', 'Arial', ...
    'FontSize', 15);
title('f', ...
    'FontName', 'Arial', ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'units','normalized',...
    'position', [-1/18,1+1/10], ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off
c