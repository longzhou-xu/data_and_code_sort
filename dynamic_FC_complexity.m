%% dynamic FC complexity
clc; clear; close all
FC_win=zeros(96,96,101);
FCdiversity_win=zeros(101,295);
FCentropy_win=zeros(101,295);
FCmean_win=zeros(101,295);
FCstd_win=zeros(101,295);
for subject=1:295
    subject
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
        'a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
        num2str(subject), '.mat']);

    Cmean(1,:) = mean(rest_HY96_ROI(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(rest_HY96_ROI(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    rest_C = (rest_HY96_ROI(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
    BOLD = rest_C; time_point=1200; ROI_number=96; win_size=200; win_number=101;
    [FC_win(:,:,:),FCdiversity_win(:,subject),FCentropy_win(:,subject),FCmean_win(:,subject),FCstd_win(:,subject),win_end] =...
        dyn_FC_FCD_FCE_FCmean_FCstd(BOLD,time_point,ROI_number,win_size,win_number);
    mkdir('FC_win');
    save(['FC_win\sub',num2str(subject),'.mat'],'FC_win');
end
save('dynamical_FC_complexity','FCdiversity_win','FCentropy_win','FCmean_win','FCstd_win','win_end');
