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
F6_a.dynamicMS = MS_win;
F6_a.dynamicFCentropy = FCentropy_win;
F6_a.FitX = min(MS_win(:)):0.01:max(MS_win(:));
F6_a.FitY = P1 .* F6_a.FitX .^ 2 + P2 .* F6_a.FitX + P3;
F6_a.FitResult = Fitresult;
for Subject=1:295
    plot(F6_a.dynamicMS(Subject,:),F6_a.dynamicFCentropy(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(F6_a.FitX,F6_a.FitY, ...
    'Color', [1,1,1], 'LineStyle', '--', 'LineWidth', 6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic FC entropy H(FC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast','FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('A','FontName','Arial','FontSize',24, ...
    'units','normalized',...
    'position',[-1/18,1+1/18],'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
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
F6_b.dynamicMS = MS_win;
F6_b.dynamicFCdiversity = FCdiversity_win;
F6_b.FitX = min(MS_win(:)):0.01:max(MS_win(:));
F6_b.FitY = P1 .* F6_b.FitX .^ 2 + P2 .* F6_b.FitX + P3;
F6_b.FitResult = Fitresult;
for Subject=1:295
    plot(F6_b.dynamicMS(Subject,:),F6_b.dynamicFCdiversity(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(F6_b.FitX, F6_b.FitY, ...
    'Color',[1,1,1],'LineStyle','--','LineWidth',6);
set(gca,'FontName','Arial','FontSize',12)
xlabel('dynamic MS <r>_n','FontName','Arial','FontSize',15)
ylabel('dynamic FC diversity D(FC)_n','FontName','Arial','FontSize',15)
legend(F1,'$fit$','Interpreter','latex','Location','southeast', ...
    'FontName','Arial','FontSize',18,'EdgeColor',[1,1,1]);
legend('boxoff')
set(gca,'color',[0.75,0.75,0.75])
set(0,'defaultfigurecolor','w') 
title('B','FontName','Arial','FontSize',24, ...
    'units','normalized',...
    'position',[-1/18,1+1/18],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
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
F6_c.dynamicMS = MS_win;
F6_c.dynamicCorrFCSC = corr_FCwin_SC;
F6_c.FitX = min(MS_win(:)):0.01:max(MS_win(:));
F6_c.FitY = P1 .* F6_c.FitX .^ 2 + P2 .* F6_c.FitX + P3;
F6_c.FitResult = Fitresult;
for Subject=1:295
    plot(F6_c.dynamicMS(Subject,:),F6_c.dynamicCorrFCSC(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(F6_c.FitX,F6_c.FitY, ...
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
title('C','FontName','Arial','FontSize',24, ...
    'units','normalized',...
    'position',[-1/18,1+1/18],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;
%% figure 6 d
AX1=subplot(2,2,4);
[Fitresult, gof_corr, output_corr] = createFit_poly2(MS_win(:), hamming_FCwin_sc(:));
P1=Fitresult.p1;
P2=Fitresult.p2;
P3=Fitresult.p3;
F6_d.dynamicMS = MS_win;
F6_d.dynamicHammingFCSC = hamming_FCwin_sc;
F6_d.FitX = min(MS_win(:)):0.01:max(MS_win(:));
F6_d.FitY = P1 .* F6_d.FitX .^ 2 + P2 .* F6_d.FitX + P3;
F6_d.FitResult = Fitresult;
for Subject=1:295
    plot(F6_d.dynamicMS(Subject,:),F6_d.dynamicHammingFCSC(Subject,:), ...
        'Marker','*','MarkerSize',2,'LineStyle','none','LineWidth',1.5);
    hold on
end
F1=plot(F6_d.FitX, F6_d.FitY, ...
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
title('D','FontName','Arial','FontSize',24, ...
    'units','normalized',...
    'position',[-1/18,1+1/18],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(AX1,'off'); AX1.LineWidth=2; hold off;
%%
save('Figure_6.mat', 'F6_a', 'F6_b', 'F6_c', 'F6_d');