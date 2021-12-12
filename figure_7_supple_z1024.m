clc; clear; close all
figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name', 'Figure 7 S4', ...
      'Position', [0.1 0.1 0.8 0.8]);
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat'])
syn(6)=[];synE(6)=[];KOP(:,6)=[];
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat'], 'subject')
subject(6)=[];
load('HCP_iq0804.mat', 'Subject')
load('HCP_iq0804.mat', 'ListSort_Unadj')
load('HCP_iq0804.mat', 'PMAT24_A_CR')
load('HCP_iq0804.mat', 'PicVocab_Unadj')
% select the subject 
for i = 1 : 294
    for j = 1 : 1206
        if subject(i) == Subject(j)
            PMAT(i, 1) = PMAT24_A_CR(j);
            LS(i, 1) = ListSort_Unadj(j);
            PV(i, 1) = PicVocab_Unadj(j);
        end
    end
end
%% figure 7 a
index_PMAT = find(PMAT > -1000);
[fit_SE_PMAT] = createFit_poly1(synE(index_PMAT,1), PMAT(index_PMAT,1));
p1 = fit_SE_PMAT.p1; 
p2 = fit_SE_PMAT.p2; 
F7_S4_a.SE = synE(index_PMAT,1);
F7_S4_a.PMAT = PMAT(index_PMAT,1);
F7_S4_a.fitX = min(synE(index_PMAT,1))-0.1:0.1:max(synE(index_PMAT,1))+0.1;
F7_S4_a.fitY = p1 .* F7_S4_a.fitX + p2;
F7_S4_a.fitResult = fit_SE_PMAT;
[R_SE_PMAT,P_SE_PMAT] = corr(F7_S4_a.SE, F7_S4_a.PMAT);

AX1 = subplot(2,3,1);
plot(F7_S4_a.SE, F7_S4_a.PMAT, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_a.fitX, F7_S4_a.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.80, 0.95, ['$R=', num2str(round(R_SE_PMAT,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize', 18);
text(0.80, 0.85, ['$p=', num2str(round(P_SE_PMAT,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize', 18);
ylim([-3 35])
xlabel('Atlas_2_4_6 SE H(r)', 'FontName', 'Arial', 'FontSize', 15);
ylabel('PMAT', 'FontName', 'Arial', 'FontSize', 15);
title('a', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units','normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off

%% figure 7 b
index_PV = find(PV > -1000);
[fit_SE_PV] = createFit_poly1(synE(index_PV,1), PV(index_PV,1));
p1 = fit_SE_PV.p1; 
p2 = fit_SE_PV.p2; 
F7_S4_b.SE = synE(index_PV,1);
F7_S4_b.PV = PV(index_PV,1);
F7_S4_b.fitX = min(synE(index_PV,1))-0.1:0.1:max(synE(index_PV,1))+0.1;
F7_S4_b.fitY = p1 .* F7_S4_b.fitX + p2;
F7_S4_b.fitResult = fit_SE_PV;
[R_SE_PV,P_SE_PV] = corr(F7_S4_b.SE, F7_S4_b.PV);

AX1 = subplot(2,3,2);
plot(F7_S4_b.SE , F7_S4_b.PV, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_b.fitX, F7_S4_b.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.80, 0.95, ['$R=', num2str(round(R_SE_PV,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
text(0.80, 0.85, ['$p=', num2str(round(P_SE_PV,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
ylim([85 165]);
xlabel('Atlas_2_4_6 SE H(r)', 'FontName', 'Arial', 'FontSize',15);
ylabel('Picture vocabulary', 'FontName', 'Arial', 'FontSize',15)
title('b', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth=2; 
hold off

%% figure 7 c
index_LS = find(LS > -1000);
[fit_SE_LS] = createFit_poly1(synE(index_LS,1), LS(index_LS,1));
p1 = fit_SE_LS.p1; 
p2 = fit_SE_LS.p2; 
F7_S4_c.SE = synE(index_LS,1);
F7_S4_c.LS = LS(index_LS,1);
F7_S4_c.fitX = min(synE(index_LS,1))-0.1:0.1:max(synE(index_LS,1))+0.1;
F7_S4_c.fitY = p1 .* F7_S4_c.fitX + p2;
F7_S4_c.fitResult = fit_SE_LS;

[R_SE_LS,P_SE_LS] = corr(F7_S4_c.SE, F7_S4_c.LS);
AX1 = subplot(2,3,3);
plot(F7_S4_c.SE, F7_S4_c.LS, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_c.fitX, F7_S4_c.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.80, 0.95, ['$R=', num2str(round(R_SE_LS,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize',18);
text(0.80, 0.85, ['$p=', num2str(round(P_SE_LS,2)), '$'], ...
    'units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontSize', 18);
ylim([60 165])
xlabel('Atlas_2_4_6 SE H(r)', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting', 'FontName', 'Arial', 'FontSize', 15)
title('c', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth = 2; 
hold off

%% figure 7 d
index_PMAT = find(PMAT > -1000);
[fit_MS_PMAT] = createFit_poly2(syn(index_PMAT,1), PMAT(index_PMAT,1));
p1 = fit_MS_PMAT.p1;
p2 = fit_MS_PMAT.p2; 
p3 = fit_MS_PMAT.p3; 
F7_S4_d.MS = syn(index_PMAT,1);
F7_S4_d.PMAT = PMAT(index_PMAT,1);
F7_S4_d.fitX = min(syn(index_PMAT,1))-0.1:0.1:max(syn(index_PMAT,1))+0.1;
F7_S4_d.fitY = p1 .* F7_S4_d.fitX .^ 2+p2 .* F7_S4_d.fitX + p3;
F7_S4_d.fitResult = fit_MS_PMAT;

AX1 = subplot(2,3,4);
plot(F7_S4_d.MS, F7_S4_d.PMAT, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_d.fitX, F7_S4_d.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 18;
lgd.Location = 'northwest';
xlim([0.1, 0.8])
ylim([-3, 35])
xlabel('Atlas_2_4_6 MS <r>', 'FontName', 'Arial', 'FontSize', 15);
ylabel('PMAT', 'FontName', 'Arial', 'FontSize', 15);
title('d', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1, 'off'); 
AX1.LineWidth = 2; 
hold off

%% figure 7 d
index_PV = find(PV > -1000);
[fit_MS_PV] = createFit_poly2(syn(index_PV,1), PV(index_PV,1));
p1 = fit_MS_PV.p1;
p2 = fit_MS_PV.p2; 
p3 = fit_MS_PV.p3; 
F7_S4_e.MS = syn(index_PV,1);
F7_S4_e.PV = PV(index_PV,1);
F7_S4_e.fitX = min(syn(index_PV,1))-0.1:0.1:max(syn(index_PV,1))+0.1;
F7_S4_e.fitY = p1 .* F7_S4_e.fitX .^ 2+p2 .* F7_S4_e.fitX + p3;
F7_S4_e.fitResult = fit_MS_PV;

AX1 = subplot(2,3,5);
plot(F7_S4_e.MS, F7_S4_e.PV, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_e.fitX, F7_S4_e.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
lgd = legend(f1, '$fit$', 'Interpreter','latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 18;
lgd.Location = 'northwest';
ylim([85 165]);
xlim([min(F7_S4_e.fitX) max(F7_S4_e.fitX)])
xlabel('Atlas_2_4_6 MS <r>', 'FontName', 'Arial', 'FontSize', 15);
ylabel('Picture vocabulary', 'FontName', 'Arial', 'FontSize', 15);
title('e', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(AX1, 'off'); 
AX1.LineWidth = 2; 
hold off
%% figure 7 f 
index_LS = find(LS > -1000);
[fit_MS_LS] = createFit_poly2(syn(index_LS,1), LS(index_LS,1));
p1 = fit_MS_LS.p1;
p2 = fit_MS_LS.p2; 
p3 = fit_MS_LS.p3; 
F7_S4_f.MS = syn(index_LS,1);
F7_S4_f.LS = LS(index_LS,1);
F7_S4_f.fitX = min(syn(index_LS,1))-0.1:0.1:max(syn(index_LS,1))+0.1;
F7_S4_f.fitY = p1 .* F7_S4_f.fitX .^ 2+p2 .* F7_S4_f.fitX + p3;
F7_S4_f.fitResult = fit_MS_LS;

AX1 = subplot(2,3,6);
plot(F7_S4_f.MS, F7_S4_f.LS, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S4_f.fitX, F7_S4_f.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize',2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 18;
lgd.Location = 'northwest';
ylim([60 165]);
xlabel('Atlas_2_4_6 MS <r>', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting', 'FontName', 'Arial', 'FontSize', 15)
title('f', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom');
grid off; 
box(AX1, 'off'); 
AX1. LineWidth = 2; 
hold off
