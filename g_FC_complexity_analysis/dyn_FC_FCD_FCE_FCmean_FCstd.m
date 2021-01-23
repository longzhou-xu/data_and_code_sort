function [FC,FCD,FCE,FCmean,FCstd,win_end] = dyn_FC_FCD_FCE_FCmean_FCstd(time_series,time_point,ROI_number,win_size,win_number)
%DYN_FCD_FCE 
%[FC_diveristy,FCE,slip] 
%  = dyn_FCD_FCE(time_series,time_point,ROI_number,win_size,win_number)
%   此处显示详细说明
win_end=floor(linspace(win_size,time_point,win_number));
%dynamcal FC
FC=zeros(ROI_number,ROI_number,win_number);
for win=1:win_number
    FC(1:ROI_number,1:ROI_number,win)=corr(time_series(win_end(win)-win_size+1:win_end(win),1:ROI_number));
end
%dynamcis FCD and FCE
for win=1:win_number
    FC_win=FC(:,:,win);
    I=0;
    FC_serise=zeros(ROI_number*(ROI_number-1)/2,1);
    for roi1=1:ROI_number
        for roi2=roi1+1:ROI_number
            I=I+1;
            FC_serise(I,1)=FC_win(roi1,roi2);
        end
    end
    %HY_96的全脑FC均值和方差
    FC_mean=mean(abs(FC_serise));                              
    FC_std=std(abs(FC_serise));                                
    %HY_96的全脑FC的熵，diversity
    b=hist(abs(FC_serise),linspace(0,1,30));
    p=b./(ROI_number*(ROI_number-1)/2);
    FC_diveristy=1;
    FD_bin_number=length(p);
    for M=1:FD_bin_number
        FC_diveristy=FC_diveristy-1./(2.*((FD_bin_number-1)./FD_bin_number)).*abs(p(M)-1./FD_bin_number);  
    end
    p(p==0)=[];
    FC_entropy=0;
    for kk=1:1:length(p)
        FC_entropy=FC_entropy-p(kk).*log2(p(kk));              
    end
    FCmean(win,1)=FC_mean;
    FCstd(win,1)=FC_std;
    FCD(win,1)=FC_diveristy;
    FCE(win,1)=FC_entropy;
end
end

