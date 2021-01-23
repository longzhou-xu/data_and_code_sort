clc; clear; close all
%% FC complexity
FC=zeros(96,96,295);
FCdiversity=zeros(295,1);
FCentropy=zeros(295,1);
FCmean=zeros(295,1);
FCstd=zeros(295,1);
for subject=1:295
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid',...
        '_Iq\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',...
        num2str(subject), '.mat']);
    Cmean(1,:) = mean(rest_HY96_ROI(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(rest_HY96_ROI(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    rest_C = (rest_HY96_ROI(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    BOLD = rest_C;
    time_point=1200;
    ROI_number=96;
    win_size=1200;
    win_number=1;
    [FC(:,:,subject),FCdiversity(subject,1),FCentropy(subject,1),...
        FCmean(subject,1),FCstd(subject,1)]=dyn_FC_FCD_FCE_FCmean_FCstd(BOLD,...
        time_point,ROI_number,win_size,win_number);
    
end
save('static_FC_complexity','FC','FCdiversity','FCentropy','FCmean','FCstd');
