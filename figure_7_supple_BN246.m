clc; clear; close all

figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name', 'Figure 7 S3', ...
      'Position', [0.1 0.1 0.8 0.8]);
  
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_BA246_RS.mat')
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat', 'subject')
load('HCP_iq0804.mat', 'Subject')
load('HCP_iq0804.mat', 'ListSort_Unadj')
load('HCP_iq0804.mat', 'PMAT24_A_CR')
load('HCP_iq0804.mat', 'PicVocab_Unadj')

% select the subject 
for i = 1 : 295
    for j = 1 : 1206
        if subject(i) == Subject(j)
            PMAT(i, 1) = PMAT24_A_CR(j);
            Listsorting(i, 1) = ListSort_Unadj(j);
            PictureVocab(i, 1) = PicVocab_Unadj(j);
        end
    end
end

%% figure 7 a
index_PMAT = find(PMAT > -1000);
[fit_SE_PMAT] = createFit_poly1(synE(index_PMAT,1), PMAT(index_PMAT,1));
p1 = fit_SE_PMAT.p1; 
p2 = fit_SE_PMAT.p2; 
F7_S3_a.SE = synE(index_PMAT,1);
F7_S3_a.PMAT = PMAT(index_PMAT,1);
F7_S3_a.fitX = min(synE(index_PMAT,1))-0.1:0.1:max(synE(index_PMAT,1))+0.1;
F7_S3_a.fitY = p1 .* F7_S3_a.fitX + p2;
F7_S3_a.fitResult = fit_SE_PMAT;
[R_SE_PMAT,P_SE_PMAT] = corr(F7_S3_a.SE, F7_S3_a.PMAT);

ax1 = subplot(2,3,1);
plot(F7_S3_a.SE, F7_S3_a.PMAT, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_a.fitX, F7_S3_a.fitY, 'Color', [1,0,0], ...
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
box(ax1,'off'); 
ax1.LineWidth = 2; 
hold off

%% figure 7 b
index_PictureVocab = find(PictureVocab > -1000);
[fit_SE_PictureVocab] = createFit_poly1(synE(index_PictureVocab,1), PictureVocab(index_PictureVocab,1));
p1 = fit_SE_PictureVocab.p1; 
p2 = fit_SE_PictureVocab.p2; 
F7_S3_b.SE = synE(index_PictureVocab,1);
F7_S3_b.PictureVocab = PictureVocab(index_PictureVocab,1);
F7_S3_b.fitX = min(synE(index_PictureVocab,1))-0.1:0.1:max(synE(index_PictureVocab,1))+0.1;
F7_S3_b.fitY = p1 .* F7_S3_b.fitX + p2;
F7_S3_b.fitResult = fit_SE_PictureVocab;

[R_SE_PictureVocab,P_SE_PictureVocab] = corr(F7_S3_b.SE, F7_S3_b.PictureVocab);
ax1 = subplot(2,3,2);
plot(F7_S3_b.SE , F7_S3_b.PictureVocab, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_b.fitX, F7_S3_b.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.80, 0.95, ['$R=', num2str(round(R_SE_PictureVocab,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
text(0.80, 0.85, ['$p=', num2str(round(P_SE_PictureVocab,2)), '$'], ...
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
box(ax1,'off'); 
ax1.LineWidth=2; 
hold off

%% figure 7 c
index_Listsorting = find(Listsorting > -1000);
[fit_SE_ListSort] = createFit_poly1(synE(index_Listsorting,1), Listsorting(index_Listsorting,1));
p1 = fit_SE_ListSort.p1; 
p2 = fit_SE_ListSort.p2; 
F7_S3_c.SE = synE(index_Listsorting,1);
F7_S3_c.ListSortingWorkMemory = Listsorting(index_Listsorting,1);
F7_S3_c.fitX = min(synE(index_Listsorting,1))-0.1:0.1:max(synE(index_Listsorting,1))+0.1;
F7_S3_c.fitY = p1 .* F7_S3_c.fitX + p2;
F7_S3_c.fitResult = fit_SE_ListSort;

[R_SE_ListSort,P_SE_ListSort] = corr(F7_S3_c.SE, F7_S3_c.ListSortingWorkMemory);
ax1 = subplot(2,3,3);
plot(F7_S3_c.SE, F7_S3_c.ListSortingWorkMemory, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_c.fitX, F7_S3_c.fitY, 'Color', [1,0,0], ...
    'LineStyle', '--', 'Linewidth', 2.5, ...
    'Marker', 'none', 'MarkerSize', 2.5);
set(gca, 'FontName', 'Arial', 'FontSize', 12)
lgd = legend(f1, '$fit$', 'Interpreter', 'latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 15;
lgd.Location = 'northwest';
text(0.80, 0.95, ['$R=', num2str(round(R_SE_ListSort,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize',18);
text(0.80, 0.85, ['$p=', num2str(round(P_SE_ListSort,2)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
    'FontName', 'Arial', 'FontSize', 18);
ylim([60 165])
xlabel('Atlas_2_4_6 SE H(r)', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting', 'FontName', 'Arial', 'FontSize', 15)
title('c', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(ax1,'off'); 
ax1.LineWidth = 2; 
hold off

%% figure 7 d
index_PMAT = find(PMAT > -1000);
[fit_MS_PMAT] = createFit_poly2(syn(index_PMAT,1), PMAT(index_PMAT,1));
p1 = fit_MS_PMAT.p1;
p2 = fit_MS_PMAT.p2; 
p3 = fit_MS_PMAT.p3; 
F7_S3_d.MS = syn(index_PMAT,1);
F7_S3_d.PMAT = PMAT(index_PMAT,1);
F7_S3_d.fitX = min(syn(index_PMAT,1))-0.1:0.1:max(syn(index_PMAT,1))+0.1;
F7_S3_d.fitY = p1 .* F7_S3_d.fitX .^ 2+p2 .* F7_S3_d.fitX + p3;
F7_S3_d.fitResult = fit_MS_PMAT;

ax1 = subplot(2,3,4);
plot(F7_S3_d.MS, F7_S3_d.PMAT, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_d.fitX, F7_S3_d.fitY, 'Color', [1,0,0], ...
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
box(ax1, 'off'); 
ax1.LineWidth = 2; 
hold off

%% figure 7 d
index_PictureVocab = find(PictureVocab > -1000);
[fit_MS_PictureVocab] = createFit_poly2(syn(index_PictureVocab,1), PictureVocab(index_PictureVocab,1));
p1 = fit_MS_PictureVocab.p1;
p2 = fit_MS_PictureVocab.p2; 
p3 = fit_MS_PictureVocab.p3; 
F7_S3_e.MS = syn(index_PictureVocab,1);
F7_S3_e.PictureVocab = PictureVocab(index_PictureVocab,1);
F7_S3_e.fitX = min(syn(index_PictureVocab,1))-0.1:0.1:max(syn(index_PictureVocab,1))+0.1;
F7_S3_e.fitY = p1 .* F7_S3_e.fitX .^ 2+p2 .* F7_S3_e.fitX + p3;
F7_S3_e.fitResult = fit_MS_PictureVocab;

ax1 = subplot(2,3,5);
plot(F7_S3_e.MS, F7_S3_e.PictureVocab, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_e.fitX, F7_S3_e.fitY, 'Color', [1,0,0], ...
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
xlim([min(F7_S3_e.fitX) max(F7_S3_e.fitX)])
xlabel('Atlas_2_4_6 MS <r>', 'FontName', 'Arial', 'FontSize', 15);
ylabel('Picture vocabulary', 'FontName', 'Arial', 'FontSize', 15);
title('e', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(ax1, 'off'); 
ax1.LineWidth = 2; 
hold off
%% figure 7 f 
index_Listsorting = find(Listsorting > -1000);
[fit_MS_Listsorting] = createFit_poly2(syn(index_Listsorting,1), Listsorting(index_Listsorting,1));
p1 = fit_MS_Listsorting.p1;
p2 = fit_MS_Listsorting.p2; 
p3 = fit_MS_Listsorting.p3; 
F7_S3_f.MS = syn(index_Listsorting,1);
F7_S3_f.Listsorting = Listsorting(index_Listsorting,1);
F7_S3_f.fitX = min(syn(index_Listsorting,1))-0.1:0.1:max(syn(index_Listsorting,1))+0.1;
F7_S3_f.fitY = p1 .* F7_S3_f.fitX .^ 2+p2 .* F7_S3_f.fitX + p3;
F7_S3_f.fitResult = fit_MS_Listsorting;

ax1 = subplot(2,3,6);
plot(F7_S3_f.MS, F7_S3_f.Listsorting, 'Color', [0.5,0.5,0.5], ...
    'LineStyle', 'none', 'Linewidth', 2, ...
    'Marker', 'o', 'MarkerSize', 5);
hold on
f1 = plot(F7_S3_f.fitX, F7_S3_f.fitY, 'Color', [1,0,0], ...
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
xlim([0.1 0.8]);
xlabel('Atlas_2_4_6 MS <r>', 'FontName', 'Arial', 'FontSize',15);
ylabel('List sorting', 'FontName', 'Arial', 'FontSize', 15)
title('f', 'FontName', 'Arial', ...
    'FontSize', 24, 'FontWeight', 'bold', ...
    'units', 'normalized', 'position', [-1/18,1+1/18], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
grid off; 
box(ax1, 'off'); 
ax1. LineWidth = 2; 
hold off
