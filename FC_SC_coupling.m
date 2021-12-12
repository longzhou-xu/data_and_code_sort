clc; clear; close all
load('SC\general_sc_series.mat')
corr_SC_FC=zeros(295,80,200);
corr_SC_FC_p=zeros(295,80,200);
hamming_SC_FC=zeros(295,80,200);
for S=1:295
    disp(S);
    load(['FC\FC_binary\sub_',num2str(S)], 'FC_binary'); %导入功能连接
    FC_binary(FC_binary>0) = 1;
    for thr=1:1:80 
        J=0;
        fc_serise=zeros(96*95/2,1);
        for roi1=1:96
            for roi2=roi1+1:96
                J=J+1;
                fc_serise(J,1)=FC_binary(roi1,roi2,thr);
            end
        end
        for s_thr=1:1:200
            [corr_SC_FC(S,thr,s_thr),corr_SC_FC_p(S,thr,s_thr)]=corr(sc_serise(:,s_thr),fc_serise);
            sum_diff=abs(sc_serise(:,s_thr)-fc_serise);
            hamming_SC_FC(S,thr,s_thr)=sum(sum_diff)./length(fc_serise);
        end
    end
end
save('SC_FC_similarity_1.mat','corr_SC_FC','corr_SC_FC_p','hamming_SC_FC')