clc;clear;close all
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_HY96_RS.mat'])
load('dynamic_MS_SE_HY96.mat')
load('MS_SE_dynamics_HY96.mat')
load('KOP_group.mat')
load('cross_interval.mat')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\HMS_MMS_LMS_subject.mat'])

figure(...
      'Color', 'w', ...
      'Units', 'Normalized', ...
      'Name', 'Figure 5 ', ...
      'Position', [0.1 0.1 0.8 0.8]);
 %% figure 5 a
[MS_SE_fit] = createFit_poly2(syn(:,1), synE(:,1));
[MS_min(1),MS_min(2)]=min(syn(:,1));
[MS_max(1),MS_max(2)]=max(syn(:,1));
P1=MS_SE_fit.p1; 
P2=MS_SE_fit.p2; 
P3=MS_SE_fit.p3;
F5_a.staticMS = syn;
F5_a.staticSE = synE;
F5_a.staticFitX = min(syn):0.01:max(syn)+0.01;
F5_a.staticFitY = P1 .* (F5_a.staticFitX) .^ 2 + P2 .* F5_a.staticFitX + P3;
F5_a.staticFitResult = MS_SE_fit;
F5_a.dynamicMS_S1 = MS_win(91,:);
F5_a.dynamicSE_S1 = SE_win(91,:);
F5_a.dynamicMS_S2 = MS_win(94,:);
F5_a.dynamicSE_S2 = SE_win(94,:);
F5_a.dynamicMS_S3 = MS_win(86,:);
F5_a.dynamicSE_S3 = SE_win(86,:);
F5_a.dynamicMS_S4 = MS_win(283,:);
F5_a.dynamicSE_S4 = SE_win(283,:);
F5_a.dynamicMS_S5 = MS_win(83,:);
F5_a.dynamicSE_S5 = SE_win(83,:);
F5_a.dynamicMS_S6 = MS_win(166,:);
F5_a.dynamicSE_S6 = SE_win(166,:);
F5_a.staticMS_S1 = syn(91,1);
F5_a.staticSE_S1 = synE(91,1);
F5_a.staticMS_S2 = syn(94,1);
F5_a.staticSE_S2 = synE(94,1);
F5_a.staticMS_S3 = syn(86,1);
F5_a.staticSE_S3 = synE(86,1);
F5_a.staticMS_S4 = syn(283,1);
F5_a.staticSE_S4 = synE(283,1);
F5_a.staticMS_S5 = syn(83,1);
F5_a.staticSE_S5 = synE(83,1);
F5_a.staticMS_S6 = syn(166,1);
F5_a.staticSE_S6 = synE(166,1);

AX1=subplot(2,3,1);
f1=plot(F5_a.staticMS, F5_a.staticSE, ...
    'color',[0.5,0.5,0.5],'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
hold on
plot(F5_a.staticFitX, F5_a.staticFitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',1);
f2=plot(F5_a.dynamicMS_S1, F5_a.dynamicSE_S1, ...
    'color',[0,0.45,0.74],'Marker','<','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
f3=plot(F5_a.dynamicMS_S2, F5_a.dynamicSE_S2, ...
    'color',[0.93,0.69,0.13],'Marker','s','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
f4=plot(F5_a.dynamicMS_S3, F5_a.dynamicSE_S3, ...
    'color',[1,0,0],'Marker','d','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
f5=plot(F5_a.dynamicMS_S4, F5_a.dynamicSE_S4, ...
    'color',[0.49,0.18,0.56],'Marker','h','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
f6=plot(F5_a.dynamicMS_S5, F5_a.dynamicSE_S5, ...
    'color',[0,0.5,0.2],'Marker','p','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
F7=plot(F5_a.dynamicMS_S6, F5_a.dynamicSE_S6, ...
    'color',[0.85,0.33,0.1],'Marker','o','MarkerSize',7,'LineStyle','none','LineWidth',1.5);
plot(F5_a.staticMS_S1,F5_a.staticSE_S1, ...
    'color',[0,0,0],'Marker','<','MarkerSize',10,'LineStyle','none','LineWidth',3);
plot(F5_a.staticMS_S2,F5_a.staticSE_S2, ...
    'color',[0,0,0],'Marker','s','MarkerSize',10,'LineStyle','none','LineWidth',3);
plot(F5_a.staticMS_S3,F5_a.staticSE_S3, ...
    'color',[0,0,0],'Marker','d','MarkerSize',10,'LineStyle','none','LineWidth',3);
plot(F5_a.staticMS_S4,F5_a.staticSE_S4, ...
    'color',[0,0,0],'Marker','h','MarkerSize',10,'LineStyle','none','LineWidth',3);
plot(F5_a.staticMS_S5,F5_a.staticSE_S5, ...
    'color',[0,0,0],'Marker','p','MarkerSize',10,'LineStyle','none','LineWidth',3);
plot(F5_a.staticMS_S6,F5_a.staticSE_S6, ...
    'color',[0,0,0],'Marker','o','MarkerSize',10,'LineStyle','none','LineWidth',3);

set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic SE H(r)_n','FontName','Arial','FontSize',15)
lgd=legend([f2,f3,f4,f5,f6,F7],...
    {'$s1$';'$s2$';'$s3$';'$s4$';'$s5$';'$s6$'},...
    'Interpreter','latex');
legend('boxoff'); 
lgd.NumColumns = 3; 
lgd.FontName = 'Arial'; 
lgd.FontSize = 18; 
lgd.Location = 'southwest';
title('A', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;
%% figure 5 b
F5_b.timepoints = 1:1200;
F5_b.KOPS5 = KOP(:,83);
F5_b.KOPS6 = KOP(:,166);
F5_b.KOPS3 = KOP(:,86);
F5_b.KOPS4 = KOP(:,283);
F5_b.KOPS1 = KOP(:,91);
F5_b.KOPS2 = KOP(:,94);
F5_b.referenceLine = 0.5*ones(1200,1);

AX1=subplot(6,3,2);
plot(F5_b.timepoints,F5_b.KOPS5, ...
    'color',[0,0.5,0.2],'Marker','p','MarkerSize',1,'LineStyle','-','LineWidth',1);
hold on
plot(F5_b.timepoints,F5_b.KOPS6, ...
    'color',[0.85,0.33,0.1],'Marker','o','MarkerSize',1,'LineStyle','-','LineWidth',1);
plot(F5_b.timepoints,F5_b.referenceLine, ...
    'color',[0,0,0],'Marker','p','MarkerSize',1,'LineStyle','-','LineWidth',1);
ylim([0,1]);
set(gca,'FontName','Arial','FontSize',12)
ylabel('r(t)','FontName','Arial','FontSize',15)
set(gca,'xtick',[],'xticklabel',[])
set(gca,'XColor','white')
title('B', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/7],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;

AX1 = subplot('Position', [0.4108 0.7014 0.2134 0.1026]);
f3=plot(F5_b.timepoints,F5_b.KOPS3, ...
    'color',[1,0,0],'Marker','d','MarkerSize',1,'LineStyle','-','LineWidth',1);
hold on
f4=plot(F5_b.timepoints,F5_b.KOPS4, ...
    'color',[0.49,0.18,0.56],'Marker','h','MarkerSize',1,'LineStyle','-','LineWidth',1);
plot(F5_b.timepoints,F5_b.referenceLine, ...
    'color',[0,0,0],'Marker','p','MarkerSize',1,'LineStyle','-','LineWidth',1);
ylim([0,1]);
set(gca,'FontName','Arial','FontSize',12)
ylabel('r(t)','FontName','Arial','FontSize',15)
set(gca,'xtick',[],'xticklabel',[])
set(gca,'XColor','white')
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;

AX1 = subplot('Position', [0.4108 0.5838 0.2134 0.1026]);
f1=plot(F5_b.timepoints,F5_b.KOPS1, ...
    'color',[0,0.45,0.74],'Marker','<','MarkerSize',1,'LineStyle','-','LineWidth',1);
hold on
f2=plot(F5_b.timepoints,F5_b.KOPS2, ...
    'color',[0.93,0.69,0.13],'Marker','s','MarkerSize',1,'LineStyle','-','LineWidth',1);
plot(F5_b.timepoints,F5_b.referenceLine, ...
    'color',[0,0,0],'Marker','p','MarkerSize',1,'LineStyle','-','LineWidth',1);
ylim([0,1]);
set(gca,'FontName','Arial','FontSize',12)
xlabel('time (vol.)','FontName','Arial','FontSize',15)
ylabel('r(t)','FontName','Arial','FontSize',15)
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;
%% figure 5 c
F5_c.histCenter=KOP_group.center_KOP_group8;
F5_c.histCount=KOP_group.count_KOP_group8;
AX1=subplot(2,3,3);
 plot(F5_c.histCenter,F5_c.histCount(:,1), ...
     'color',[1,0,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 hold on
 plot(F5_c.histCenter,F5_c.histCount(:,2), ...
     'color',[0,1,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,3), ...
     'color',[0,0,1],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,4), ...
     'color',[0.5,0.5,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,5), ...
     'color',[0.5,0.5,1],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,6), ...
     'color',[0,0.5,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,7), ...
     'color',[1,0.5,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 plot(F5_c.histCenter,F5_c.histCount(:,8), ...
     'color',[0.5,0,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
 
set(gca,'FontName','Arial','FontSize',12);
lgd = legend('$<0.30$', '$[0.30,0.35)$', '$[0.35,0.40)$', ...
    '$[0.40,0.45)$', '$[0.45,0.50)$', '$[0.50,0.55)$', ...
    '$[0.55,0.60)$', '$>0.60$',...
    'Interpreter','latex');
legend('boxoff'); 
lgd.NumColumns = 2; 
lgd.FontName = 'Arial'; 
lgd.FontSize = 12; 
lgd.Position = [0.7594 0.8672 0.1636 0.0953];
xlabel('r(t)', 'FontName', 'Arial', 'FontSize', 15);
ylabel('PDF', 'FontName', 'Arial', 'FontSize', 15);
title('C', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; 
box(AX1, 'off'); 
axis tight; 
AX1.LineWidth = 2; 
hold off;

%% figure 5 d
clear count center
F5_d.histcount = cross_interval.count_cf_interval_group8;
F5_d.histcenter = cross_interval.center_cf_interval_group8;
J=1;
AX1=subplot(2,3,4);
f1=loglog(F5_d.histcenter{1,J},F5_d.histcount{1,J}, ...
    'color',[1,0,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
hold on
f2=loglog(F5_d.histcenter{2,J},F5_d.histcount{2,J}, ...
    'color',[0,1,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f3=loglog(F5_d.histcenter{3,J},F5_d.histcount{3,J}, ...
    'color',[0,0,1],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f4=loglog(F5_d.histcenter{4,J},F5_d.histcount{4,J}, ...
    'color',[0.5,0.5,0],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f5=loglog(F5_d.histcenter{5,J},F5_d.histcount{5,J}, ...
    'color',[0.5,0.5,1],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f6=loglog(F5_d.histcenter{6,J},F5_d.histcount{6,J}, ...
    'color',[0,0.5,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f7=loglog(F5_d.histcenter{7,J},F5_d.histcount{7,J}, ...
    'color',[1,0.5,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
f8=loglog(F5_d.histcenter{8,J},F5_d.histcount{8,J}, ...
    'color',[0.5,0,0.5],'Marker','none','MarkerSize',1,'LineStyle','-','LineWidth',1);
xlim([1,200])
set(gca,'FontName','Arial','FontSize',12);
lgd=legend('$<0.30$','$[0.30,0.35)$','$[0.35,0.40)$','$[0.40,0.45)$','$[0.45,0.50)$','$[0.50,0.55)$','$[0.55,0.60)$','$>0.60$',...
    'Interpreter','latex');
legend('boxoff');lgd.NumColumns = 1;lgd.FontName = 'Arial';lgd.FontSize = 15;lgd.Location = 'southwest';
xlabel('T (vol.)','FontName','Arial','FontSize',15);
ylabel('PDF','FontName','Arial','FontSize',15);
title('D', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off;
box(AX1,'off');
AX1.LineWidth=2;
hold off;
%% figure 5 e
AX1 = subplot(4,3,8);
[F5_e.histcount(:,1),F5_e.histcenter(:,1)] = ...
    hist(reshape(dynamics.move_syn_win,100*295,1),100);
F5_e.histcount(:,1) = F5_e.histcount ./ sum(F5_e.histcount);
bar(F5_e.histcenter(:,1), F5_e.histcount(:,1));
set(gca,'FontName','Arial','FontSize',12)
xlabel('\Delta<r>_n','FontName','Arial','FontSize',15)
ylabel('PDF','FontName','Arial','FontSize',15)
title('E', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2;
%% figure 5 f
AX1=subplot(4,3,9);
[F5_f.histcount(:,1),F5_f.histcenter(:,1)] = hist(reshape(dynamics.move_synE_win,100*295,1),100);
F5_f.histcount(:,1) = F5_f.histcount ./ sum(F5_f.histcount);
bar(F5_f.histcenter(:,1), F5_f.histcount(:,1));
set(gca,'FontName','Arial','FontSize',12)
xlabel('\DeltaH(r)_n','FontName','Arial','FontSize',15)
ylabel('PDF','FontName','Arial','FontSize',15)
title('F', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');grid off
box(AX1,'off');
AX1.LineWidth=2;
%% figure 5 g
AX1=subplot(4,3,11);
[MS_MSspeed_fit] = createFit_poly2(syn(:,1), mean(dynamics.speed_syn_win,2)./(10*0.72));
P1=MS_MSspeed_fit.p1;
P2=MS_MSspeed_fit.p2;
P3=MS_MSspeed_fit.p3;

F5_g.staticMS = syn(:,1);
F5_g.dynamicMSspeed = mean(dynamics.speed_syn_win,2)./(10*0.72);
F5_g.fitX = min(syn):0.01:max(syn)+0.01;
F5_g.fitY = P1 .* (F5_g.fitX) .^ 2 + P2 .* F5_g.fitX + P3;
F5_g.fitResult = MS_MSspeed_fit;

plot(F5_g.staticMS,F5_g.dynamicMSspeed, ...
    'color',[0.5,0.5,0.5],'Marker','*','MarkerSize',3,'LineStyle','none','LineWidth',1.5);
hold on
f1=plot(F5_g.fitX,F5_g.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2.5);
lgd=legend(f1,'$fit$','Interpreter','latex');
legend('boxoff');
lgd.NumColumns = 1;
lgd.FontName = 'Arial';
lgd.FontSize = 12;
lgd.Location = 'northwest';
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('<|\Delta<r>_n|>/\Deltan','Interpreter','tex','FontName','Arial','FontSize',15)
title('G', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off
spss{3,1}=syn(:,1);
spss{3,2}=mean(dynamics.speed_syn_win,2)./(10*0.72);
%%
AX1=subplot(4,3,12);
[MS_SEspeed_fit] = createFit_poly2(syn(:,1), mean(dynamics.speed_synE_win,2)./(10*0.72));
P1=MS_SEspeed_fit.p1;
P2=MS_SEspeed_fit.p2;
P3=MS_SEspeed_fit.p3;
clear x y
x=min(syn):0.01:max(syn)+0.01;
y=P1.*(x).^2+P2.*x+P3;

F5_h.staticMS = syn(:,1);
F5_h.dynamicSEspeed = mean(dynamics.speed_synE_win,2)./(10*0.72);
F5_h.fitX = min(syn):0.01:max(syn)+0.01;
F5_h.fitY = P1 .* (F5_h.fitX) .^ 2 + P2 .* F5_h.fitX + P3;
F5_h.fitResult = MS_SEspeed_fit;

plot(F5_h.staticMS,F5_h.dynamicSEspeed, ...
    'color',[0.5,0.5,0.5],'Marker','*','MarkerSize',3,'LineStyle','none','LineWidth',1.5);
hold on
f1=plot(F5_h.fitX,F5_h.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2.5);
set(gca,'FontName','Arial','FontSize',12)
lgd=legend(f1,'$fit$','Interpreter','latex');
legend('boxoff');lgd.NumColumns = 1;lgd.FontName = 'Arial';lgd.FontSize = 12;lgd.Location = 'northwest';
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('<|\DeltaH(r)_n|>/\Deltan','Interpreter','tex','FontName','Arial','FontSize',15)
title('H', 'FontName','Arial',...
    'FontSize',24, ...
    'units','normalized', 'position',[-1/18,1+1/18],...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2;
spss{4,1}=syn(:,1);
spss{4,2}=mean(dynamics.speed_synE_win,2)./(10*0.72);

save('Figure_5.mat','F5_a','F5_b','F5_c','F5_d','F5_e','F5_f','F5_g','F5_h')
