clc; clear; close all

load('SC_FC_similarity.mat')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_HY96_RS.mat']);
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\HMS_MMS_LMS_subject.mat']);
figure
thr=70;
s_thr=40;
[fitresult_corr] = createFit_poly2(syn, corr_SC_FC(:,thr,s_thr));
Spss(:,1)=syn;
Spss(:,2)=corr_SC_FC(:,thr,s_thr);
p1=fitresult_corr.p1;
p2=fitresult_corr.p2;
p3=fitresult_corr.p3;

figure_4_a.MS = syn;
figure_4_a.PearsonFC_SC = corr_SC_FC(:,thr,s_thr);
figure_4_a.fitX = min(syn(:,1)):0.01:max(syn(:,1));
figure_4_a.fitY = p1 .* figure_4_a.fitX .^ 2 + p2 .* figure_4_a.fitX + p3;
figure_4_a.fitResult = fitresult_corr;

ax1=subplot(2,2,1);
plot(figure_4_a.MS, figure_4_a.PearsonFC_SC, ...
    'color',[0.5,0.5,0.5],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
hold on
f1=plot(figure_4_a.fitX, figure_4_a.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12);
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC) ','FontName','Arial','FontSize',15)
text(0.7,0.96,['$\rho_{FC}=',num2str(thr/100),'$'],'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
text(0.7,0.86,['$THR_{SC}=',num2str(s_thr),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
legend(f1,'$fit$','Interpreter','latex','Location','northwest','FontName','Arial','FontSize',17,'EdgeColor',[1,1,1]);
legend('boxoff')
title('a','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units',...
    'normalized','position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; ax1.LineWidth=2; box(ax1,'off'); hold off

ax1=subplot(2,2,2);
[fitresult_Hamming, gof_H, output_H] = createFit_poly2(syn, hamming_SC_FC(:,thr,s_thr));
Spss(:,3)=hamming_SC_FC(:,thr,s_thr);
p1=fitresult_Hamming.p1;
p2=fitresult_Hamming.p2;
p3=fitresult_Hamming.p3;

figure_4_b.MS = syn;
figure_4_b.hamming_SC_FC = hamming_SC_FC(:,thr,s_thr);
figure_4_b.fitX = min(syn(:,1)):0.01:max(syn(:,1));
figure_4_b.fitY = p1 .* figure_4_b.fitX .^ 2 + p2 .* figure_4_b.fitX + p3;
figure_4_b.fitResult = fitresult_Hamming;

plot(figure_4_b.MS, figure_4_b.hamming_SC_FC, ...
    'color',[0.5,0.5,0.5],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
hold on
f1=plot(figure_4_b.fitX,figure_4_b.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15);
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15);
text(0.7,0.96,['$\rho_{FC}=',num2str(thr/100),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00], ...
        'FontName','Arial','FontSize',17,'EdgeColor','none');
text(0.7,0.86,['$THR_{SC}=',num2str(s_thr),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00], ...
        'FontName','Arial','FontSize',17,'EdgeColor','none');
legend(f1,'$fit$','Interpreter','latex','Location','southwest','FontName','Arial','FontSize',17,'EdgeColor',[1,1,1]);
legend('boxoff')
title('b','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized','position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; ax1.LineWidth=2; box(ax1,'off'); hold off;

ax1=subplot(2,2,3);
HMS_corr=corr_SC_FC(suppercri_subject_number(:,3),:,s_thr);
MMS_corr=corr_SC_FC(cri_subject_number(:,3),:,s_thr);
LMS_corr=corr_SC_FC(subcri_subject_number(:,3),:,s_thr);
for thr=1:80
    [h23(thr,1),p23(thr,1), ~, stats23{thr,1}]=ttest2(HMS_corr(:,thr),MMS_corr(:,thr));
    [h12(thr,1),p12(thr,1), ~, stats12{thr,1}]=ttest2(LMS_corr(:,thr),MMS_corr(:,thr));
end

figure_4_c.FCdensityRange = 0.01:0.01:0.8;
figure_4_c.meanHMScorrFCSC = mean(HMS_corr,1);
figure_4_c.meanMMScorrFCSC = mean(MMS_corr,1);
figure_4_c.meanLMScorrFCSC = mean(LMS_corr,1);
figure_4_c.ttest2_HMS_MMS_stats = stats23;
figure_4_c.ttest2_LMS_MMS_stats = stats12;
figure_4_c.ttest2_HMS_MMS_p = p23;
figure_4_c.ttest2_LMS_MMS_p = p12;

f1=plot(figure_4_c.FCdensityRange, figure_4_c.meanHMScorrFCSC, ...
    'Color',[0.85,0.33,0.1],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
hold on
f2=plot(figure_4_c.FCdensityRange,figure_4_c.meanMMScorrFCSC, ...
    'Color',[0.0,0.45,0.74],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
f3=plot(figure_4_c.FCdensityRange,figure_4_c.meanLMScorrFCSC, ...
    'Color',[0.93,0.69,0.13],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
[p23_005_index(:,1)]=find(p23<0.05);
[p12_005_index(:,1)]=find(p12<0.05);
plot(p23_005_index./100,figure_4_c.meanHMScorrFCSC(p23_005_index)+0.002,'Color',[0.49,0.18,0.56],'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(p12_005_index./100,figure_4_c.meanLMScorrFCSC(p12_005_index)+0.002,'Color',[0.47,0.67,0.19],'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);
set(gca,'FontName','Arial','FontSize',12)
xlabel('FC density \rho_F_C','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
text(0.7,0.96,['$THR_{SC}=',num2str(s_thr),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
legend([f1,f2,f3],'$HMS$','$MMS$','$LMS$','Interpreter','latex','Location','southwest',...
    'FontName','Arial','FontSize',17,'EdgeColor',[1,1,1]);
legend('boxoff')
title('c','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized','position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off
ax1.LineWidth=2;
box(ax1,'off')  
hold off

ax1=subplot(2,2,4);
HMS_hamming=hamming_SC_FC(suppercri_subject_number(:,3),:,s_thr);
MMS_hamming=hamming_SC_FC(cri_subject_number(:,3),:,s_thr);
LMS_hamming=hamming_SC_FC(subcri_subject_number(:,3),:,s_thr);
clear h23 h23 h12 p12 stats23 stats12
for thr=1:80
    [h23(thr,1),p23(thr,1), ~, stats23{thr,1}]=ttest2(HMS_hamming(:,thr),MMS_hamming(:,thr));
    [h12(thr,1),p12(thr,1), ~, stats12{thr,1}]=ttest2(LMS_hamming(:,thr),MMS_hamming(:,thr));
end

figure_4_d.FCdensityRange = 0.01:0.01:0.8;
figure_4_d.meanHMScorrFCSC = mean(HMS_hamming,1);
figure_4_d.meanMMScorrFCSC = mean(MMS_hamming,1);
figure_4_d.meanLMScorrFCSC = mean(LMS_hamming,1);
figure_4_d.ttest2_HMS_MMS_stats = stats23;
figure_4_d.ttest2_LMS_MMS_stats = stats12;
figure_4_d.ttest2_HMS_MMS_p = p23;
figure_4_d.ttest2_LMS_MMS_p = p12;

f1=plot(figure_4_d.FCdensityRange, figure_4_d.meanHMScorrFCSC, ...
    'Color',[0.85,0.33,0.1],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
hold on
f2=plot(figure_4_d.FCdensityRange, figure_4_d.meanMMScorrFCSC, ...
    'Color',[0.0,0.45,0.74],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
f3=plot(figure_4_d.FCdensityRange, figure_4_d.meanLMScorrFCSC, ...
    'Color',[0.93,0.69,0.13],'LineStyle','-','LineWidth',2,'Marker','none','MarkerSize',3);
clear p23_005_index p12_005_index
[p23_005_index(:,1)]=find(p23<0.05);
[p12_005_index(:,1)]=find(p12<0.05);
plot(p23_005_index./100, figure_4_d.meanHMScorrFCSC(p23_005_index)+0.003, ...
    'Color',[0.49,0.18,0.56],'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(p12_005_index./100, figure_4_d.meanLMScorrFCSC(p12_005_index)+0.003, ...
    'Color',[0.47,0.67,0.19],'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);
set(gca,'FontName','Arial','FontSize',12)
xlabel('FC density \rho_F_C','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
text(0.7,0.96,['$THR_{SC}=',num2str(s_thr),'$'],'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
legend('$HMS$','$MMS$','$LMS$','Interpreter','latex','Location','northwest', ...
    'FontName','Arial','FontSize',17,'EdgeColor',[1,1,1]);
legend('boxoff')
title('d','FontName','Arial','FontSize',24, 'FontWeight', 'bold', ...
    'units','normalized','position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; ax1.LineWidth=2; box(ax1,'off'); hold off

save('figure_4_a.mat', 'figure_4_a');
save('figure_4_b.mat', 'figure_4_b');
save('figure_4_c.mat', 'figure_4_c');
save('figure_4_d.mat', 'figure_4_d');