clc;clear;close all
% Calculate the MS and SE of HY96
for S = 1:295
    S
    for JJ = 2:100
        JJ
       % read the data
       load(['D:\supplyment_dataANDcode\a_extract_the_signals_from_HCP_files\ROI_signals\HY_96\sub',num2str(S),'.mat']);
       signal=rest_HY96_ROI(:,1:96);
       % Z-score normalization
       Cmean(1,:) = mean(signal(1:1200,:),1);
       Cmeanmatrix = repmat(Cmean,1200,1);
       Cstd = std(signal(1:1200,:),1);
       Cstdmatrix = repmat(Cstd,1200,1);
       signal_z_score = (signal(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
       %calculate the MS and SE
       [syn(S,JJ),synE(S,JJ),KOP(:,S)] = syn_synEntropy(signal_z_score(1:1200,:),1200,96,JJ);
    end
end
save('MS_SE_HY96_RS_S1.mat','syn','synE');
%%
clc;clear;close all
load('MS_SE_HY96_RS_S1.mat', 'syn')
load('MS_SE_HY96_RS_S1.mat', 'synE')
figure
subplot(3,3,1)
plot(syn(:,5), synE(:,5),'o');
subplot(3,3,2)
plot(syn(:,15), synE(:,15),'o');
subplot(3,3,3)
plot(syn(:,25), synE(:,25),'o');
subplot(3,3,4)
plot(syn(:,35), synE(:,35),'o');
subplot(3,3,5)
plot(syn(:,45), synE(:,45),'o');
subplot(3,3,6)
plot(syn(:,55), synE(:,55),'o');
subplot(3,3,7)
plot(syn(:,65), synE(:,65),'o');
subplot(3,3,8)
plot(syn(:,75), synE(:,75),'o');
subplot(3,3,9)
plot(syn(:,100), synE(:,100),'o');

figure
% for II =1:100
%     plot(syn(:,II), synE(:,II),'o');
%     pause(1)
% end