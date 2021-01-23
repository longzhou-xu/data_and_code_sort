clc;clear;close all;
for S=1:295
    S
    load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\',...
        'a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S)]);
    % normalize the signal
    Cmean(1,:) = mean(rest_HY96_ROI(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(rest_HY96_ROI(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    Rest_C = (rest_HY96_ROI(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    Surrogate_data = phaseran_lz(Rest_C(1:1200,:),500);
    mkdir('surrogate_data');
    save(['surrogate_data\surrogate_sub_',num2str(S),'.mat'],'Surrogate_data');
end