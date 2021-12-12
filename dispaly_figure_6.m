clc
clear 
close all
load('dynamical_FC_complexity.mat', 'FCdiversity_win')
load('dynamical_FC_complexity.mat', 'FCentropy_win')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_HY96_RS.mat'], 'syn')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\d_MS_SE_dynamic\dynamic_MS_SE_HY96.mat'], 'MS_win')
FCentropy_win=FCentropy_win';
FCdiversity_win=FCdiversity_win';
%% figure 6 a
AX1=subplot(2,2,1);
Spss(:,1)=MS_win(:); 
Spss(:,2)=FCentropy_win(:);
[Fitresult] = createFit_poly2(MS_win(:), FCentropy_win(:));
P1=Fitresult.p1; 
P2=Fitresult.p2; 
P3=Fitresult.p3;

Figure_6_a.dynamicMS = MS_win;
Figure_6_a.dynamicFCentropy = FCentropy_win;
Figure_6_a.FitX = min(MS_win(:)):0.01:max(MS_win(:));
Figure_6_a.FitY = P1 .* Figure_6_a.FitX .^ 2 + P2 .* Figure_6_a.FitX + P3;
Figure_6_a.FitResult = Fitresult;
save('Figure_6_a.mat','Figure_6_a');

for Subject=1:295
    plot(Figure_6_a.dynamicMS(Subject,:),Figure_6_a.dynamicFCentropy(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(Figure_6_a.FitX,Figure_6_a.FitY, ...
    'Color', [1,1,1], 'LineStyle', '--', 'LineWidth', 6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic FC entropy H(FC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast','FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('a','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; 
box(AX1,'off'); 
AX1.LineWidth=2; 
hold off;
%% figure 6 b
AX1=subplot(2,2,2);
[Fitresult] = createFit_poly2(MS_win(:), FCdiversity_win(:));
P1=Fitresult.p1; 
P2=Fitresult.p2; 
P3=Fitresult.p3;

Figure_6_b.dynamicMS = MS_win;
Figure_6_b.dynamicFCdiversity = FCdiversity_win;
Figure_6_b.FitX = min(MS_win(:)):0.01:max(MS_win(:));
Figure_6_b.FitY = P1 .* Figure_6_b.FitX .^ 2 + P2 .* Figure_6_b.FitX + P3;
Figure_6_b.FitResult = Fitresult;
save('Figure_6_b.mat','Figure_6_b');
Spss(:,3)=FCdiversity_win(:);
for Subject=1:295
    plot(Figure_6_b.dynamicMS(Subject,:),Figure_6_b.dynamicFCdiversity(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(Figure_6_b.FitX, Figure_6_b.FitY, ...
    'Color',[1,1,1],'LineStyle','--','LineWidth',6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic FC diversity D(FC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast', ...
    'FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('b','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;

%% figure 6 c and d
thr=7;
s_thr=4;
for Subject=1:295
    load(['dynamic_FC_SC_couple\sub',num2str(Subject),'.mat'])
    corr_FCwin_SC(Subject,:)=corr_SC_FCwin(:,thr,s_thr);
    hamming_FCwin_sc(Subject,:)=hamming_SC_FCwin(:,thr,s_thr);
end

%% figure 6 c
AX1=subplot(2,2,3);

[Fitresult] = createFit_poly2(MS_win(:), corr_FCwin_SC(:));
P1=Fitresult.p1;
P2=Fitresult.p2;
P3=Fitresult.p3;
Figure_6_c.dynamicMS = MS_win;
Figure_6_c.dynamicCorrFCSC = corr_FCwin_SC;
Figure_6_c.FitX = min(MS_win(:)):0.01:max(MS_win(:));
Figure_6_c.FitY = P1 .* Figure_6_c.FitX .^ 2 + P2 .* Figure_6_c.FitX + P3;
Figure_6_c.FitResult = Fitresult;
save('Figure_6_c.mat','Figure_6_c');
Spss(:,4)=corr_FCwin_SC(:);

for Subject=1:295
    plot(Figure_6_c.dynamicMS(Subject,:),Figure_6_c.dynamicCorrFCSC(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(Figure_6_c.FitX,Figure_6_c.FitY, ...
    'Color',[1,1,1],'LineStyle','--','LineWidth',6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic R(FC-SC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast', ...
    'FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
text(0.75,0.95,['$\rho_{FC}=',num2str(thr/10),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
text(0.75,0.85,['$THR_{SC}=',num2str(s_thr*10),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('c','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;

%% figure 6 d
AX1=subplot(2,2,4);
[Fitresult, gof_corr, output_corr] = createFit_poly2(MS_win(:), hamming_FCwin_sc(:));
P1=Fitresult.p1;
P2=Fitresult.p2;
P3=Fitresult.p3;
Figure_6_d.dynamicMS = MS_win;
Figure_6_d.dynamicHammingFCSC = hamming_FCwin_sc;
Figure_6_d.FitX = min(MS_win(:)):0.01:max(MS_win(:));
Figure_6_d.FitY = P1 .* Figure_6_d.FitX .^ 2 + P2 .* Figure_6_d.FitX + P3;
Figure_6_d.FitResult = Fitresult;
save('Figure_6_d.mat','Figure_6_d');
Spss(:,4)=corr_FCwin_SC(:);

for Subject=1:295
    plot(Figure_6_d.dynamicMS(Subject,:),Figure_6_d.dynamicHammingFCSC(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(Figure_6_d.FitX, Figure_6_d.FitY, ...
    'Color',[1,1,1],'LineStyle','--','LineWidth',6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic HD(FC-SC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast', ...
    'FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
text(0.75,0.95,['$\rho_{FC}=',num2str(thr/10),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
text(0.75,0.85,['$THR_{SC}=',num2str(s_thr*10),'$'],'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('d','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/18],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;
