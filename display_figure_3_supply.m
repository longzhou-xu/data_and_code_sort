clc; clear; close all;

load('flexibility_supple.mat', 'CNE_subject')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_HY96_RS.mat']);
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_Z1024_RS.mat'], 'subject')


CNE=CNE_subject(:,10,11);%选择展示thr=0.3，win number=60的结果，其他作为补充材料

index_1 = find(subject == 101107);
index_2 = find(subject == 134728);
index_include = 1:295;
index_include(index_include == index_2)=[];
index_include(index_include == index_1)=[];

%% supplement figure S3 
figure
[Syn_sort,Syn_sort_index]=sort(syn(:,1));
Range=[0,14,51,106,155,202,254,275,295];
for I=1:length(Range)-1
    Syn_mean(I,1)=mean(Syn_sort(Range(I)+1:Range(I+1)));
    for Thr=0.12:0.02:0.48  %(0.12:0.02:0.48)
        Bin_num=60; %(40:2:60)
        clear Cne_60
        X=int8((Thr-0.1)./0.02);
        Y=(Bin_num-38)./2;
        Cne_60(:,1)=CNE_subject(Syn_sort_index(Range(I)+1:Range(I+1)),X,Y);
        Cne_60_mean(I,X)=mean(Cne_60);
    end
end
for I=1:length(Range)-1
    for Bin_num=40:2:60
        Thr=0.3;
        clear Cne_03
        X=int8((Thr-0.1)./0.02);
        Y=(Bin_num-38)./2;
        Cne_03(:,1)=CNE_subject(Syn_sort_index(Range(I)+1:Range(I+1)),X,Y);
        Cne_03_mean(I,Y)=mean(Cne_03);
    end
end
AX1=subplot(1,2,1);
plot(Syn_mean,Cne_60_mean(:,1),'Marker','o','MarkerSize',6,'LineStyle','-','LineWidth',2);
hold on
plot(Syn_mean,Cne_60_mean(:,4),'Marker','+','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_60_mean(:,7),'Marker','*','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_60_mean(:,10),'Marker','.','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_60_mean(:,13),'Marker','x','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_60_mean(:,16),'Marker','s','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_60_mean(:,19),'Marker','d','MarkerSize',6,'LineStyle','-','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('flexibility CNE','FontName','Arial','FontSize',15)
legend('$THR_{FC}=0.12$','$THR_{FC}=0.18$','$THR_{FC}=0.24$','$THR_{FC}=0.30$','$THR_{FC}=0.36$',...
    '$THR_{FC}=0.42$','$THR_{FC}=0.48$','Interpreter','latex',...
    'Location','northwest','FontName','Arial','FontSize',20,'EdgeColor',[1,1,1]);
legend('boxoff')
title('a','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized','position',[-1/18,1+1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

ax4=subplot(1,2,2);
plot(Syn_mean,Cne_03_mean(:,1),'Marker','o','MarkerSize',6,'LineStyle','-','LineWidth',2);
hold on
plot(Syn_mean,Cne_03_mean(:,3),'Marker','+','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_03_mean(:,5),'Marker','*','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_03_mean(:,7),'Marker','.','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_03_mean(:,9),'Marker','s','MarkerSize',6,'LineStyle','-','LineWidth',2);
plot(Syn_mean,Cne_03_mean(:,11),'Marker','d','MarkerSize',6,'LineStyle','-','LineWidth',2);
set(gca,'FontName','Arial','FontSize',15)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('flexibility CNE','FontName','Arial','FontSize',15)
legend('$n_{win}=40$','$n_{win}=44$','$n_{win}=48$','$n_{win}=52$','$n_{win}=56$','$n_{win}=60$',...
    'Interpreter','latex','Location','northwest','FontName','Arial','FontSize',20,'EdgeColor',[1,1,1]);
legend('boxoff')
title('b','FontName','Arial','FontSize',24, 'FontWeight', 'bold', 'units','normalized','position',[-1/18,1+1/30],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; ax4.LineWidth=2; box(ax4,'off'); hold off;
