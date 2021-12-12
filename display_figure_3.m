clc
clear 
close all
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
    'b_MS_SE_STATIC\MS_SE_HY96_RS.mat'], 'syn')
load('flexibility.mat', 'CNE_subject')
load('static_FC_complexity.mat', 'FCdiversity')
load('static_FC_complexity.mat', 'FCentropy')
load('static_FC_complexity.mat', 'FCmean')

CNE=CNE_subject;
index_1 = 6;
index_2 = 173;
index_include = 1:295;
index_include(index_include == index_2)=[];
index_include(index_include == index_1)=[];
spss1(:,1)=syn(:,1);
spss1(:,2)=FCmean(:,1);
spss1(:,3)=FCentropy(:,1);
spss1(:,4)=FCdiversity(:,1);
spss1(:,5)=CNE(:,1);
spss2(:,1)=syn(index_include,1);
spss2(:,2)=FCmean(index_include,1);
spss2(:,3)=FCentropy(index_include,1);
spss2(:,4)=FCdiversity(index_include,1);
spss2(:,5)=CNE(index_include,1);
%% Figure 3 
AX1=subplot(1,3,1);
[fit_MS_FCE] = createFit_poly2(syn(index_include,1), FCentropy(index_include,1));
p1=fit_MS_FCE.p1;
p2=fit_MS_FCE.p2;
p3=fit_MS_FCE.p3;
clear x y
F3_a.MSincludeCorr = syn(index_include,1);
F3_a.FCentropyIncludeCorr = FCentropy(index_include);
F3_a.MSoutlier_1 = syn(index_1,1);
F3_a.FCentropyOutlier_1 = FCentropy(index_1);
F3_a.MSoutlier_2 = syn(index_2,1);
F3_a.FCentropyOutlier_2 = FCentropy(index_2);
F3_a.fitX = min(syn):0.01:max(syn)+0.01;
F3_a.fitY = p1 .* F3_a.fitX .^ 2+p2 .* F3_a.fitX + p3;
F3_a.fitResult = fit_MS_FCE;

plot(F3_a.MSincludeCorr, F3_a.FCentropyIncludeCorr, ...
    'color',[0.5,0.5,0.5],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
hold on
f1=plot(F3_a.MSoutlier_1, F3_a.FCentropyOutlier_1, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
plot(F3_a.MSoutlier_2, F3_a.FCentropyOutlier_2, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
f2=plot(F3_a.fitX, F3_a.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('FC entropy H(FC)','FontName','Arial','FontSize',15)
legend([f1,f2],'$outliers$','$fit$','Interpreter','latex', ...
    'Location','southeast','FontName','Arial','FontSize',15,'EdgeColor',[1,1,1]);
legend('boxoff')
title('a','FontName','Arial','FontSize',24, 'FontWeight', 'bold', ...
    'units','normalized','position',[-1/18,1+1/80],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

AX1=subplot(1,3,2);
[fit_MS_FCD] = createFit_poly2(syn(index_include,1), FCdiversity(index_include,1));
p1=fit_MS_FCD.p1;
p2=fit_MS_FCD.p2;
p3=fit_MS_FCD.p3;

F3_b.MSincludeCorr = syn(index_include,1);
F3_b.FCdiversityIncludeCorr = FCdiversity(index_include);
F3_b.MSoutlier_1 = syn(index_1,1);
F3_b.FCdiversityOutlier_1 = FCdiversity(index_1);
F3_b.MSoutlier_2 = syn(index_2,1);
F3_b.FCdiversityOutlier_2 = FCdiversity(index_2);
F3_b.fitX = min(syn):0.01:max(syn)+0.01;
F3_b.fitY = p1 .* F3_b.fitX .^ 2+p2 .* F3_b.fitX + p3;
F3_b.fitResut = fit_MS_FCD;

plot(F3_b.MSincludeCorr, F3_b.FCdiversityIncludeCorr, ...
    'color',[0.5,0.5,0.5],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
hold on
f1=plot(F3_b.MSoutlier_1, F3_b.FCdiversityOutlier_1, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
plot(F3_b.MSoutlier_2, F3_b.FCdiversityOutlier_2, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
f2=plot(F3_b.fitX, F3_b.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12);
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('FC diversity D(FC)','FontName','Arial','FontSize',15)
legend([f1,f2],'$outliers$','$fit$','Interpreter','latex', ...
    'Location','southeast','FontName','Arial','FontSize',15,'EdgeColor',[1,1,1]);
legend('boxoff')
title('b','FontName','Arial','FontSize',24, 'FontWeight', 'bold', ...
    'units','normalized','position',[-1/18,1+1/80],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2; box(AX1,'off'); hold off

AX1=subplot(1,3,3);
connect_number_entropy1_include=CNE(index_include,1);
[fit_MS_CNE] = createFit_poly2(syn(index_include,1), CNE(index_include,1));
p1=fit_MS_CNE.p1;
p2=fit_MS_CNE.p2;
p3=fit_MS_CNE.p3;

F3_c.MSincludeCorr = syn(index_include,1);
F3_c.CNEincludeCorr = CNE(index_include,1);
F3_c.MSoutlier_1 = syn(index_1,1);
F3_c.CNEoutlier_1 = CNE(index_1);
F3_c.MSoutlier_2 = syn(index_2,1);
F3_c.CNEoutlier_2 = CNE(index_2);
F3_c.fitX = min(syn):0.01:max(syn)+0.01;
F3_c.fitY = p1 .* F3_c.fitX .^ 2+p2 .* F3_c.fitX + p3;
F3_c.fitResut = fit_MS_CNE;

plot(F3_c.MSincludeCorr,F3_c.CNEincludeCorr, ...
    'color',[0.5,0.5,0.5],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
hold on
f1=plot(F3_c.MSoutlier_1,F3_c.CNEoutlier_1, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
plot(F3_c.MSoutlier_2,F3_c.CNEoutlier_2, ...
    'color',[1,0,0],'Marker','o','MarkerSize',6,'LineStyle','none','LineWidth',2);
f2=plot(F3_c.fitX,F3_c.fitY, ...
    'Color',[1,0,0],'LineStyle','--','LineWidth',2);
set(gca,'FontName','Arial','FontSize',12)
xlabel('MS <r>','FontName','Arial','FontSize',15)
ylabel('flexibility CNE','FontName','Arial','FontSize',15)
legend([f1,f2],'$outliers$','$fit$','Interpreter','latex', ...
    'Location','southeast','FontName','Arial','FontSize',15,'EdgeColor',[1,1,1]);
title('c','FontName','Arial','FontSize',24, 'FontWeight', 'bold', ...
    'units','normalized','position',[-1/18,1+1/80],'HorizontalAlignment','right','VerticalAlignment','bottom');
grid off; AX1.LineWidth=2;box(AX1,'off');hold off

save('figure_3_a.mat','F3_a');
save('figure_3_b.mat','F3_b');
save('figure_3_c.mat','F3_c');