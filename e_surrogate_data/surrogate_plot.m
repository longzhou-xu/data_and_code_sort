clc;clear;close all
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub1.mat'])
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'e_MS_SE_surrogate\surrogate_noremain_FC\surrogate_data\',...
    'surrogate_sub_1.mat'])
load('static_MS_SE_HY96_RS_surrogate.mat')
Surrogate_lz_data=Surrogate_data;
syn_surr_lz = syn;
synE_surr_lz = synE;

Signal=Surrogate_lz_data(:,:,1);
SZ=size(Signal);
Bmean(1,:) = mean(Signal(1:SZ(1),:),1);
Bmeanmatrix = repmat(Bmean,SZ(1),1);
Bstd = std(Signal(1:SZ(1),:),1);
Bstdmatrix = repmat(Bstd,SZ(1),1);
rest_B = (Signal(1:SZ(1),:)-Bmeanmatrix)./Bstdmatrix;

Signal=rest_HY96_ROI;
SZ=size(Signal);
Cmean(1,:) = mean(Signal(1:SZ(1),:),1);
Cmeanmatrix = repmat(Cmean,SZ(1),1);
Cstd = std(Signal(1:SZ(1),:),1);
Cstdmatrix = repmat(Cstd,SZ(1),1);
rest_C = (Signal(1:SZ(1),:)-Cmeanmatrix)./Cstdmatrix;
figure
set(gca,'color',[1,1,1])
set(0,'defaultfigurecolor','w') 
ax1=subplot(1,2,1);
plot(1:100,rest_C(1:100,1),'color',[1,0,0],'Marker','none','MarkerSize',6,'LineStyle','-','LineWidth',1.5);
hold on
plot(1:100,rest_B(1:100,1),'color',[0,0,1],'Marker','none','MarkerSize',6,'LineStyle','-','LineWidth',1.5);
set(gca,'FontName','Arial','FontSize',18);
xlabel('time (vol.)','FontName','Arial','FontSize',18)
ylabel('BOLD','FontName','Arial','FontSize',18)
lgd=legend('$surrogate$','$real$','Interpreter','latex');
legend('boxoff'); lgd.NumColumns = 1; lgd.FontName = 'Arial'; lgd.FontSize = 15; lgd.Location = 'northwest';
title('a','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/100],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(ax1,'off'); ax1.LineWidth=2; hold off;

ax1=subplot(1,2,2);
for I=1:500
   plot(syn_surr_lz(:,I),synE_surr_lz(:,I),'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',1.5);
   hold on
end
set(gca,'FontName','Arial','FontSize',18);
xlabel('MS <r>','FontName','Arial','FontSize',18)
ylabel('SE H(r)','FontName','Arial','FontSize',18)
title('b','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized',...
    'position',[-1/18,1+1/100],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; box(ax1,'off'); ax1.LineWidth=2; hold off;
