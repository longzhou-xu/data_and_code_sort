clc; clear; close all
%---------------load the MS and SE-------------------------------------
load('D:\Criticality_cognitive\b_MS_SE_STATIC\MS_SE_HY96_RS.mat','syn','synE');
load('D:\criticality_cognitive\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat', 'subject');
Subject_MS_SE = subject;clear subject;
%---------------load the FC complexity------------------------------------
load('D:\Criticality_cognitive\g_FC_complexity_analysis\static_FC_complexity.mat', 'FCdiversity');
load('D:\Criticality_cognitive\g_FC_complexity_analysis\static_FC_complexity.mat', 'FCentropy');
load('D:\Criticality_cognitive\g_FC_complexity_analysis\flexibility.mat', 'CNE_subject');
%---------------load the HCP cognitive scores-----------------------------
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
if sum(Subject_age_edu == Subject_cognitive)~=1206
    error('the subject ID is wrong');
end
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
               PMAT(XX, 8) = FCentropy(II);
               PMAT(XX, 9) = FCdiversity(II);
               PMAT(XX, 10) = CNE_subject(II);
               %-----------------LS-----------------
               LS(XX, 1) = II;
               LS(XX, 2) = ListSort_Unadj(JJ);
               LS(XX, 3) = age(JJ);
               LS(XX, 4) = edu(JJ);
               LS(XX, 5) = syn(II);
               LS(XX, 6) = synE(II);
               LS(XX, 7) = Subject_MS_SE(II);
               LS(XX, 8) = FCentropy(II);
               LS(XX, 9) = FCdiversity(II);
               LS(XX, 10) = CNE_subject(II);
               %-----------------PV----------------
               PV(XX, 1) = II;
               PV(XX, 2) = PicVocab_Unadj(JJ);
               PV(XX, 3) = age(JJ);
               PV(XX, 4) = edu(JJ);
               PV(XX, 5) = syn(II);
               PV(XX, 6) = synE(II);
               PV(XX, 7) = Subject_MS_SE(II);
               PV(XX, 8) = FCentropy(II);
               PV(XX, 9) = FCdiversity(II);
               PV(XX, 10) = CNE_subject(II);
            end
        end
    end
end
%% -----correlation test using the adjusted cognitive score-----------
%------------------------------for PMAT-----------------------------------
X = []; X = [ones(length(PMAT(:,1)),1),PMAT(:,3),PMAT(:,4)];
[b_PMAT] = regress(PMAT(:,2),X);
Residual_PMAT = PMAT(:,2) - X * b_PMAT; 
%-----------------for LS--------------------------------------
X = []; X = [ones(length(LS(:,1)),1),LS(:,3),LS(:,4)];
[b_LS] = regress(LS(:,2),X);
Residual_LS = LS(:,2) - X * b_LS; 
%--------------------for PV-------------------------------------
X = []; X = [ones(length(PV(:,1)),1),PV(:,3),PV(:,4)];
[b_PV] = regress(PV(:,2),X);
Residual_PV = PV(:,2) - X * b_PV; 
%%
figure('Color', 'w', 'Units', 'Normalized', 'Name', 'Figure 7 S1', 'Position', [0 0 1 1]);
%-------------------------------------------------------------------------
plot_linear(PMAT(:,8),PMAT(:,2), [2,3,1],'FC entropy H(FC)','PMAT','a');
[R,P] = corr(PMAT(:,8),PMAT(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([0,30])
title('a', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(PMAT(:,9),PMAT(:,2), [2,3,2],'FC diversity D(FC)','PMAT','a');
[R,P] = corr(PMAT(:,9),PMAT(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([0,30])
title('b', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(PMAT(:,10),PMAT(:,2), [2,3,3],'FC flexibility CNE','PMAT','b');
[R,P] = corr(PMAT(:,10),PMAT(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([0,30])
title('c', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%-------------------------------------------------------------------------
plot_linear(LS(:,8),LS(:,2), [2,3,4],'FC entropy H(FC)','List sorting','a');
[R,P] = corr(LS(:,8),LS(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([75,160])
title('d', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(LS(:,9),LS(:,2), [2,3,5],'FC diversity D(FC)','List sorting','a');
[R,P] = corr(LS(:,9),LS(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([75,160])
title('e', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(LS(:,10),LS(:,2), [2,3,6],'FC flexibility CNE','List sorting','b');
[R,P] = corr(LS(:,10),LS(:,2),'type','Pearson');
text(0.7, 0.98, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.88, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
ylim([75,160])
title('f', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%%
figure('Color', 'w', 'Units', 'Normalized', 'Name', 'Figure 7 SM', 'Position', [0 0 1 1]);
%-------------------------------------------------------------------------
plot_linear(PMAT(:,8),Residual_PMAT, [2,3,1],'FC entropy H(FC)','adjusted PMAT','a');
[R,P] = corr(PMAT(:,8),Residual_PMAT,'type','Pearson');
text(0.7, 1, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.9, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-15,15]);
title('a', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(PMAT(:,9),Residual_PMAT, [2,3,2],'FC diversity D(FC)','adjusted PMAT','a');
[R,P] = corr(PMAT(:,9),Residual_PMAT,'type','Pearson');
text(0.7, 1, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.9, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-15,15]);
title('b', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(PMAT(:,10),Residual_PMAT, [2,3,3],'FC flexibility CNE','adjusted PMAT','b');
[R,P] = corr(PMAT(:,10),Residual_PMAT,'type','Pearson');
text(0.7, 1, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.9, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-15,15]);
title('c', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%-------------------------------------------------------------------------
plot_linear(LS(:,8),Residual_LS, [2,3,4],'FC entropy H(FC)','adjusted List sorting','a');
[R,P] = corr(LS(:,8),Residual_LS,'type','Pearson');
text(0.7, 1, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.90, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-40,50]);
title('d', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(LS(:,9),Residual_LS, [2,3,5],'FC diversity D(FC)','adjusted List sorting','a');
[R,P] = corr(LS(:,9),Residual_LS,'type','Pearson');
text(0.7, 1.00, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.7, 0.90, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-40,50]);
title('e', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');
%--------------------------------------------------------------------------
plot_linear(LS(:,10),Residual_LS, [2,3,6],'FC flexibility CNE','adjusted List sorting','b');
[R,P] = corr(LS(:,10),Residual_LS,'type','Pearson');
text(0.75, 1, ['$R=', num2str(round(R,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
text(0.75, 0.9, ['$p=', num2str(round(P,3)), '$'], ...
    'units', 'normalized', 'Interpreter', 'latex', ...
        'FontName', 'Arial', 'FontSize', 18);
    ylim([-40,50]);
title('f', 'FontName', 'Arial', 'FontSize', 24,'FontWeight','bold');


