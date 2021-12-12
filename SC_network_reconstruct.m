clc; clear; close all
load('SC\general_sc.mat')
for s_thr=1:1:200
    GSC_1=GSC;
    GSC_1(GSC_1<s_thr)=0;
    GSC_1(GSC_1>=s_thr)=1;
    spars_gsc(s_thr,1)=sum(sum(GSC_1))/(96*96);
    J=0;
    for roi1=1:96
        for roi2=roi1+1:96
            J=J+1;
            sc_serise(J,s_thr)=GSC_1(roi1,roi2);
        end
    end
end
save('SC\general_sc_series','sc_serise','spars_gsc');