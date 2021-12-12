clc; clear; close all

figure('Color', 'w', 'Units', 'Normalized', ...
      'Name', 'Figure 7 S2', 'Position', [0.1 0.1 0.8 0.8]);
  
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_HY96_RS.mat'])
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_Z1024_RS.mat'], 'subject')
load('HCP_iq0804.mat', 'Subject')
load('HCP_iq0804.mat', 'ListSort_AgeAdj')
load('HCP_iq0804.mat', 'PicVocab_AgeAdj')

% select the corresponding subject 
for i = 1 : 295
    for j = 1 : 1206
        if subject(i) == Subject(j)
            LS(i, 1) = ListSort_AgeAdj(j);
            PV(i, 1) = PicVocab_AgeAdj(j);
        end
    end
end

%% Figure 7 S2 a
index_PV = find(PV > -1000);
[fit_SE_PV] = createFit_poly1(synE(index_PV,1), PV(index_PV,1));
p1 = fit_SE_PV.p1; 
p2 = fit_SE_PV.p2; 
F7_S2_a.SE = synE(index_PV,1);
F7_S2_a.PV = PV(index_PV,1);
F7_S2_a.fitX = min(synE(index_PV,1))-0.1:0.1:max(synE(index_PV,1))+0.1;
F7_S2_a.fitY = p1 .* F7_S2_a.fitX + p2;
F7_S2_a.fitResult = fit_SE_PV;

[R_SE_PV,P_SE_PV] = corr(F7_S2_a.SE, F7_S2_a.PV);
AX1 = subplot(2,2,1);
plot(F7_S2_a.SE , F7_S2_a.PV, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S2_a.fitX, F7_S2_a.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.1, 0.8, ['$R=', num2str(round(R_SE_PV,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
text(0.1, 0.7, ['$p=', num2str(round(P_SE_PV,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
xlabel('SE H(r)', 'FontName', 'Arial', 'FontSize',15);
ylabel('Picture vocabulary Age-adjusted', 'FontName', 'Arial', 'FontSize',15)
title('a', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth=2; 
hold off

%% Figure 7 S2 b
index_LS = find(LS > -1000);
[fit_SE_LS] = createFit_poly1(synE(index_LS,1), LS(index_LS,1));
p1 = fit_SE_LS.p1; 
p2 = fit_SE_LS.p2; 
F7_S2_b.SE = synE(index_LS,1);
F7_S2_b.LS = LS(index_LS,1);
F7_S2_b.fitX = min(synE(index_LS,1))-0.1:0.1:max(synE(index_LS,1))+0.1;
F7_S2_b.fitY = p1 .* F7_S2_b.fitX + p2;
F7_S2_b.fitResult = fit_SE_LS;

[R_SE_LS,P_SE_LS] = corr(F7_S2_b.SE, F7_S2_b.LS);
AX1 = subplot(2,2,2);
plot(F7_S2_b.SE, F7_S2_b.LS, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S2_b.fitX, F7_S2_b.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', ...
    'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.10, 0.80, ['$R=', num2str(round(R_SE_LS,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
text(0.10, 0.7, ['$p=', num2str(round(P_SE_LS,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize', 18);
xlabel('SE H(r)', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting Age-adjusted', 'FontName', 'Arial', 'FontSize', 15)
title('b', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off
%% Figure 7 S2 c
index_PV = find(PV > -1000);
[fit_MS_PV] = createFit_poly2(syn(index_PV,1), PV(index_PV,1));
p1 = fit_MS_PV.p1;
p2 = fit_MS_PV.p2; 
p3 = fit_MS_PV.p3; 
F7_S2_c.MS = syn(index_PV,1);
F7_S2_c.PV = PV(index_PV,1);
F7_S2_c.fitX = min(syn(index_PV,1))-0.1:0.1:max(syn(index_PV,1))+0.1;
F7_S2_c.fitY = p1 .* F7_S2_c.fitX .^ 2+p2 .* F7_S2_c.fitX + p3;
F7_S2_c.fitResult = fit_MS_PV;

AX1 = subplot(2,2,3);
plot(F7_S2_c.MS, F7_S2_c.PV, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S2_c.fitX, F7_S2_c.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
lgd = legend(f1, '$fit$', 'Interpreter','latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 18;
lgd.Location = 'northwest';
xlim([min(F7_S2_c.fitX) max(F7_S2_c.fitX)])
xlabel('MS <r>', 'FontName', 'Arial', 'FontSize', 15);
ylabel('Picture vocabulary Age-adjusted', 'FontName', 'Arial', 'FontSize', 15);
title('c', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1, 'off'); 
AX1.LineWidth = 2; 
hold off
%% Figure 7 S2 d
index_LS = find(LS > -1000);
[fit_MS_LS] = createFit_poly2(syn(index_LS,1), LS(index_LS,1));
p1 = fit_MS_LS.p1;
p2 = fit_MS_LS.p2; 
p3 = fit_MS_LS.p3; 
F7_S2_d.MS = syn(index_LS,1);
F7_S2_d.LS = LS(index_LS,1);
F7_S2_d.fitX = min(syn(index_LS,1))-0.1:0.1:max(syn(index_LS,1))+0.1;
F7_S2_d.fitY = p1 .* F7_S2_d.fitX .^ 2+p2 .* F7_S2_d.fitX + p3;
F7_S2_d.fitResult = fit_MS_LS;

AX1 = subplot(2,2,4);
plot(F7_S2_d.MS, F7_S2_d.LS, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S2_d.fitX, F7_S2_d.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize',2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 18;
lgd.Location = 'northwest';
xlabel('MS <r>', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting Age-adjusted', 'FontName', 'Arial', 'FontSize', 15)
title('d', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized','position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1, 'off'); 
AX1. LineWidth = 2; 
hold off
