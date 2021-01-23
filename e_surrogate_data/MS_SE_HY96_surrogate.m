clc;clear;close all
% Calculate the MS and SE of HY96 surrogate data
for S = 1:295
    S
    % read the data
    load(['surrogate_data\surrogate_sub_',num2str(S),'.mat'])
    for II = 1:500
        II
        signal=Surrogate_data(:,1:96,II);
        % Z-score normalization
        Cmean(1,:) = mean(signal(1:1199,:),1);
        Cmeanmatrix = repmat(Cmean,1199,1);
        Cstd = std(signal(1:1199,:),1);
        Cstdmatrix = repmat(Cstd,1199,1);
        signal_z_score = (signal(1:1199,:)-Cmeanmatrix)./Cstdmatrix;
        %calculate the MS and SE
        [syn(S,II),synE(S,II),KOP(:,1)] = syn_synEntropy(signal_z_score(1:1199,:),1199,96,30);
        T = 1199;
        Win_size = 200;
        Slide = 10;
        Entropy_bin = 30;
        [MS_win(S,II,:),SE_win(S,II,:)] = syn_synE_win(KOP,T,Win_size,Slide,Entropy_bin);
    end
    save('static_MS_SE_HY96_RS_surrogate.mat','syn','synE');
    save('dynamic_MS_SE_HY96_surrogate.mat', 'MS_win','SE_win');
end
