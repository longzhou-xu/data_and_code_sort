clc; clear; close all
%---------------load the MS and SE-------------------------------------
load('D:\Criticality_cognitive\b_MS_SE_STATIC\MS_SE_HY96_RS.mat');
load('D:\criticality_cognitive\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat', 'subject');
Subject_MS_SE = subject;clear Subject;
%---------------load the HCP cognitive scores---------------------------
load('HCP_iq0804.mat', 'Subject');
Subject_cognitive = Subject;clear Subject;
load('HCP_iq0804.mat', 'ListSort_Unadj');
load('HCP_iq0804.mat', 'PMAT24_A_CR');
load('HCP_iq0804.mat', 'PicVocab_Unadj');
%--------------load the HCP age and edu scores----------------------------
load('subject_age_edu.mat', 'Age_in_Yrs');
age = Age_in_Yrs(2:end);clear Age_in_Yrs;
load('subject_age_edu.mat', 'SSAGA_Educ');
edu = SSAGA_Educ(2:end);clear SSAGA_Educ;
load('subject_age_edu.mat', 'Subject');
Subject_age_edu = Subject(2:end);clear Subject;
sum(Subject_age_edu == Subject_cognitive);
%--------------check and select the subject----------------------------
%-----------demand:(1)same ID; (2)score exits------------------------------
XX =0;
for II = 1 : 295
    for JJ = 1 : 1206
        if Subject_MS_SE(II) == Subject_cognitive(JJ)
            if PMAT24_A_CR(JJ) > 0 && age(JJ) > 0 && edu(JJ) > 0 ...
                    && ListSort_Unadj(JJ) > 0 && PicVocab_Unadj(JJ) > 0
                XX = XX+1;
               %------------PMAT------------------------
               PMAT(XX, 1) = II;
               PMAT(XX, 2) = PMAT24_A_CR(JJ);
               PMAT(XX, 3) = age(JJ);
               PMAT(XX, 4) = edu(JJ);
               PMAT(XX, 5) = syn(II);
               PMAT(XX, 6) = synE(II);
               PMAT(XX, 7) = Subject_MS_SE(II);
               %-----------------LS-----------------
               LS(XX, 1) = II;
               LS(XX, 2) = ListSort_Unadj(JJ);
               LS(XX, 3) = age(JJ);
               LS(XX, 4) = edu(JJ);
               LS(XX, 5) = syn(II);
               LS(XX, 6) = synE(II);
               LS(XX, 7) = Subject_MS_SE(II);
               %-----------------PV----------------
               PV(XX, 1) = II;
               PV(XX, 2) = PicVocab_Unadj(JJ);
               PV(XX, 3) = age(JJ);
               PV(XX, 4) = edu(JJ);
               PV(XX, 5) = syn(II);
               PV(XX, 6) = synE(II);
               PV(XX, 7) = Subject_MS_SE(II);
            end
        end
    end
end
%% -----correlation test using the adjusted cognitive score-----------
%------------------------------for PMAT-----------------------------------
X = []; X = [ones(length(PMAT(:,1)),1),PMAT(:,3),PMAT(:,4)];
[b_PMAT] = regress(PMAT(:,2),X);
Residual_PMAT = PMAT(:,2) - X * b_PMAT; 
[R_SE_PMAT_adjusted,P_SE_PMAT_adjuested] = corr(Residual_PMAT,synE(PMAT(:,1)));
[R_MS_PMAT_adjusted,P_MS_PMAT_adjuested] = corr(Residual_PMAT,syn(PMAT(:,1)));
%-----------------for LS--------------------------------------
X = []; X = [ones(length(LS(:,1)),1),LS(:,3),LS(:,4)];
[b_LS] = regress(LS(:,2),X);
Residual_LS = LS(:,2) - X * b_LS; 
[R_SE_LS_adjusted,P_SE_LS_adjuested] = corr(Residual_LS,synE(LS(:,1)));
[R_MS_LS_adjusted,P_MS_LS_adjuested] = corr(Residual_LS,syn(LS(:,1)));
%--------------------for PV-------------------------------------
X = []; X = [ones(length(PV(:,1)),1),PV(:,3),PV(:,4)];
[b_PV] = regress(PV(:,2),X);
Residual_PV = PV(:,2) - X * b_PV; 
[R_SE_PV_adjusted,P_SE_PV_adjuested] = corr(Residual_PV,synE(PV(:,1)));
[R_MS_PV_adjusted,P_MS_PV_adjuested] = corr(Residual_PV,syn(PV(:,1)));
%% ---------------plot--------------------------------
figure('Color', 'w', 'Units', 'Normalized', 'Name', 'Figure 7 adjusted', 'Position', [0 0 1 1]);

plot_linear(synE(PMAT(:,1)),Residual_PMAT(:,1), [2,3,1],'SE H(r)','adjusted PMAT','a');
[R,P] = corr(synE(PMAT(:,1)),Residual_PMAT(:,1),'type','Pearson');
text(0.7, 0.95, ['$R=', num2str(round(R,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.85, ['$p=', num2str(round(P,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
xlim([min(synE(PMAT(:,1))) - 0.05,max(synE(PMAT(:,1))) + 0.05]);
ylim([-20,15])
title('a', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');

plot_linear(synE(PV(:,1)),Residual_PV(:,1),[2,3,2],'SE H(r)','adjusted Picture vocabulary','b');
[R,P] = corr(synE(PV(:,1)),Residual_PV(:,1),'type','Pearson');
text(0.7, 0.95, ['$R=', num2str(round(R,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.85, ['$p=', num2str(round(P,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
xlim([min(synE(PMAT(:,1))) - 0.05,max(synE(PMAT(:,1))) + 0.05]);
ylim([-25,45])
title('b', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');

plot_linear(synE(LS(:,1)),Residual_LS(:,1), [2,3,3],'SE H(r)','adjusted List sorting','c');
[R,P] = corr(synE(LS(:,1)),Residual_LS(:,1),'type','Pearson');
text(0.7, 0.95, ['$R=', num2str(round(R,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.85, ['$p=', num2str(round(P,3)), '$'], 'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
xlim([min(synE(PMAT(:,1))) - 0.05,max(synE(PMAT(:,1))) + 0.05]);
ylim([-45,50])
title('c', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');

plot_quadratic(syn(PMAT(:,1)),Residual_PMAT(:,1), [2,3,4],'MS <r>','adjusted PMAT','d');
xlim([min(syn(PMAT(:,1))) - 0.05,max(syn(PMAT(:,1))) + 0.05]);
ylim([-20,15])
title('d', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');

plot_quadratic(syn(PV(:,1)),Residual_PV(:,1), [2,3,5],'MS <r>','adjusted Picture vocabulary','e');
xlim([min(syn(PMAT(:,1))) - 0.05,max(syn(PMAT(:,1))) + 0.05]);
ylim([-25,45])
title('e', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');

plot_quadratic(syn(LS(:,1)),Residual_LS(:,1), [2,3,6],'MS <r>','adjusted List sorting','f');
xlim([min(syn(PMAT(:,1))) - 0.05,max(syn(PMAT(:,1))) + 0.05]);
ylim([-45,50])
title('f', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%% bar plot
figure
bin_number = 8;
A = [2,3,1];
XLABEL = 'SE H(r)';YLABEL = 'adjusted PMAT';TITLE = 'a';
barplot_SE(PMAT(:,6),Residual_PMAT,bin_number,A,XLABEL,YLABEL,TITLE);
A = [2,3,2];
XLABEL = 'SE H(r)';YLABEL = 'adjusted Picture vocabulary';TITLE = 'b';
barplot_SE(PV(:,6),Residual_PV,bin_number,A,XLABEL,YLABEL,TITLE);
A = [2,3,3];
XLABEL = 'SE H(r)';YLABEL = 'adjusted List sorting';TITLE = 'c';
barplot_SE(LS(:,6),Residual_LS,bin_number,A,XLABEL,YLABEL,TITLE);
A = [2,3,4];
XLABEL = 'MS <r>';YLABEL = 'adjusted PMAT';TITLE = 'd';
barplot_MS(PMAT(:,5),Residual_PMAT,bin_number,A,XLABEL,YLABEL,TITLE);
A = [2,3,5];
XLABEL = 'MS <r>';YLABEL = 'adjusted Picture vocabulary';TITLE = 'e';
barplot_MS(PV(:,5),Residual_PV,bin_number,A,XLABEL,YLABEL,TITLE);
A = [2,3,6];
XLABEL = 'MS <r>';YLABEL = 'adjusted List sorting';TITLE = 'f';
barplot_MS(LS(:,5),Residual_LS,bin_number,A,XLABEL,YLABEL,TITLE);