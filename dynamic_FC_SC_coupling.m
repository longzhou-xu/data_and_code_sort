clc; clear; close all
load('SC\general_sc_series.mat','sc_serise','spars_gsc');
for S=1:295
    disp(S);
    corr_SC_FCwin=zeros(101,8,15);
    corr_SC_FCwin_p=zeros(101,8,15);
    hamming_SC_FCwin=zeros(101,8,15);
    for win=1:101
        load(['FC_win_binary\sub',num2str(S),'\win',num2str(win),'\binary_FC.mat'], 'FC_binary'); 
        FC_binary(FC_binary>0) = 1;
        for thr=10:10:80 
            J=0;
            fc_serise=zeros(96*95/2,1);
            for roi1=1:96
                for roi2=roi1+1:96
                    J=J+1;
                    fc_serise(J,1)=FC_binary(roi1,roi2,thr);
                end
            end
            for s_thr=10:10:150
                [corr_SC_FCwin(win,thr/10,s_thr/10),corr_SC_FCwin_p(win,thr/10,s_thr/10)]=corr(sc_serise(:,s_thr),fc_serise);
                sum_diff=abs(sc_serise(:,s_thr)-fc_serise);
                hamming_SC_FCwin(win,thr/10,s_thr/10)=sum(sum_diff)./length(fc_serise);
            end
        end
    end
    mkdir('dynamic_FC_SC_couple')
    save(['dynamic_FC_SC_couple\sub',num2str(S),'.mat'],'corr_SC_FCwin','corr_SC_FCwin_p','hamming_SC_FCwin')
end

