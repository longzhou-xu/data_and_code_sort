clc;clear;close all
load(['D:\criticality_cognitive\project_synchrony_',...
    'avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_HY96_RS.mat'],...
    'syn');
load('SC_FC_similarity.mat')
[MS_sort,MS_sort_index]=sort(syn);
%[min, 0.3, 0.35, 0.4, 0.45, 0.5,0.55,0.6,max]
Range=[0,12,48,99,145,192,244,265,284];
MS_mean=zeros(length(Range)-1,1);
for I=1:length(Range)-1
    MS_mean(I,1)=mean(MS_sort(Range(I)+1:Range(I+1)));
end
for I=1:length(Range)-1
    for Thr=1:80
        for S_thr=1:200
            corr_sc_fc_group(I,Thr,S_thr)=mean(corr_SC_FC(MS_sort_index(Range(I)+1:Range(I+1)),Thr,S_thr));
            hanming_dis_group(I,Thr,S_thr)=mean(hamming_SC_FC(MS_sort_index(Range(I)+1:Range(I+1)),Thr,S_thr));
        end
    end
end
%% Figure 4 supplyment 1
figure
S_thr=40;
AX1=subplot(1,4,1);
for Thr=45:5:80
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',17,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$\rho_{FC}=',num2str(Thr/100),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.7])
text(0.07,0.97,['$THR_{SC}=',num2str(S_thr),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('	A','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

AX1=subplot(1,4,2);
for Thr=45:5:80
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$\rho_{FC}=',num2str(Thr/100),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.7])
text(0.07,0.97,['$THR_{SC}=',num2str(S_thr),'$'],'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('B','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

AX1=subplot(1,4,3);
for Thr=5:8:40
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$\rho_{FC}=',num2str(Thr/100),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.7])
text(0.07,0.97,['$THR_{SC}=',num2str(S_thr),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('C','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

AX1=subplot(1,4,4);
for Thr=5:8:40
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$\rho_{FC}=',num2str(Thr/100),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.7])
text(0.07,0.97,['$THR_{SC}=',num2str(S_thr),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('D','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off');  hold off

%% ²¹³äÍ¼Æ¬¶þ
figure
Thr=70;
AX1=subplot(1,4,1);
for S_thr=10:30:200
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('A','FontName','Arial','FontSize',24, ...
     'units','normalized','position',[-1/12,1-1/30],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

Thr=55;
AX1=subplot(1,4,2);
for S_thr=10:35:200
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('B','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1+1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

Thr=40;
AX1=subplot(1,4,3);
for S_thr=10:30:200
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('C','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

Thr=30;
AX1=subplot(1,4,4);
for S_thr=10:30:200
    plot(MS_mean,corr_sc_fc_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,corr_sc_fc_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('R(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('D','FontName','Arial','FontSize',24, 'units','normalized','position',[-1/12,1-1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off
%% Figure 4 supplyment 3 
figure
Thr=70;
AX1=subplot(1,4,1);
for S_thr=10:30:200
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',17,'EdgeColor','none');
title('A','FontName','Arial','FontSize',24, ...
    'units','normalized','position',[-1/12,1-1/30],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

Thr=55;
AX1=subplot(1,4,2);
for S_thr=10:40:200
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized',...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
title('B','FontName','Arial','FontSize',24, ...
    'units','normalized','position',[-1/12,1-1/30],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off');  hold off

Thr=40;
AX1=subplot(1,4,3);
for S_thr=10:30:200
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
ylim([0.36,0.5])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
title('C','FontName','Arial','FontSize',24, ...
    'units','normalized','position',[-1/12,1-1/30],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off');  hold off

Thr=30;
AX1=subplot(1,4,4);
for S_thr=10:30:200
    plot(MS_mean,hanming_dis_group(:,Thr,S_thr),'Marker','.','MarkerSize',15,'LineStyle','-','LineWidth',2);
    text(MS_mean(8)+0.01,hanming_dis_group(8,Thr,S_thr),['$THR_{SC}=',num2str(S_thr),'$'],...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',14,'EdgeColor','none');
    hold on
end
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('HD(FC-SC)','FontName','Arial','FontSize',15)
xlim([0.25,0.75])
ylim([0.32,0.51])
text(0.1,0.97,['$\rho_{FC}=',num2str(Thr/100),'$'], 'units', 'normalized', ...
        'Interpreter','latex','Color',[0.00,0.00,0.00],'FontName','Arial','FontSize',15,'EdgeColor','none');
title('D','FontName','Arial','FontSize',24,...
    'units','normalized','position',[-1/12,1-1/30],...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off');  hold off
