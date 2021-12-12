%%
clc; clear; close all
% load the KOP form static mean sychrony and synchrony entropy
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_HY96_RS.mat')
% calculate the dynamic mean synchrony and synchrony entropy
T = 1200;
Win_size = 200;
Slide = 10;
Entropy_bin = 30;
[MS_win,SE_win] = syn_synE_win(KOP,T,Win_size,Slide,Entropy_bin);
save('dynamic_MS_SE_HY96.mat', 'MS_win','SE_win');

%%
clear
% load the data 
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_HY96_RS.mat')
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\HMS_MMS_LMS_subject.mat')
% group the subject
[MS_sort,MS_sort_index] = sort(syn(:,1));
Range8 = [0,14,51,106,155,202,254,275,295];
MS_sort(Range8(2:9))
KOP_group8=cell(length(Range8)-1,1);
% select the KOP belong to the corresponding group
for I=1:length(Range8)-1 % traverse the groups
    for J=Range8(I)+1:Range8(I+1) % traverse the subject belong to the group
        KOP_group8{I,1}=[KOP_group8{I,1},KOP(:,MS_sort_index(J))];
    end
end
% hist the KOP of each group
for I=1:length(KOP_group8)
    sz=size(KOP_group8{I,1});
    [count_KOP_group8(:,I),center_KOP_group8]=hist(reshape(KOP_group8{I,1},sz(1)*sz(2),1),linspace(0,1,30));
    count_KOP_group8(:,I)=count_KOP_group8(:,I)./(sz(1)*sz(2));
end
% select the KOP belong to MMS/LMS/HMS
for I=1:20
    HS_KOP(:,I)=KOP(:,suppercri_subject_number(I,3));
    MS_KOP(:,I)=KOP(:,cri_subject_number(I,3));
    LS_KOP(:,I)=KOP(:,subcri_subject_number(I,3));
end
% save
KOP_group.KOP_group8=KOP_group8;
KOP_group.count_KOP_group8=count_KOP_group8;
KOP_group.center_KOP_group8=center_KOP_group8;
KOP_group.HS_KOP=HS_KOP;
KOP_group.MS_KOP=MS_KOP;
KOP_group.LS_KOP=LS_KOP;
save('KOP_group.mat','KOP_group');

% the across event of KOP
Thr = 0.5; % the threshold 
CF = zeros(295,1200,length(Thr));
% define the cross event
for Sub=1:295
    for I=1:length(Thr)
        for T=1:1199
            if KOP(T,Sub)<=Thr(I) && KOP(T+1,Sub)>Thr(I)
               CF(Sub,T,I)=1;
            end
        end
    end
end
% defind the cross time
for Sub=1:295
    for I=1:length(Thr)
        cf_time{Sub,I} = find(CF(Sub,:,I) == 1);
    end
end
% group the cross interval in 8 group
cf_interval_group8=cell(length(Range8)-1,length(Thr));
for I=1:length(Range8)-1
    for K=1:length(Thr)
        for J=Range8(I)+1:Range8(I+1)
            cf_interval_group8{I,K}=[cf_interval_group8{I,K},diff(cf_time{MS_sort_index(J),K})];
        end
    end
end
% group the cross interval in MMS/LMS/HMS
cf_interval_group3=cell(3,9);
for I=1:length(Thr)
    for J=1:20
        cf_interval_group3{1,I}=[cf_interval_group3{1,I},diff(cf_time{subcri_subject_number(J,3),I})];
    end
end
for I=1:length(Thr)
    for J=1:20
        cf_interval_group3{2,I}=[cf_interval_group3{2,I},diff(cf_time{cri_subject_number(J,3),I})];
    end
end
for I=1:length(Thr)
    for J=1:20
        cf_interval_group3{3,I}=[cf_interval_group3{3,I},diff(cf_time{suppercri_subject_number(J,3),I})];
    end
end
% hist the cross interval
for I=1:8
    for J=1:1
        [count_cf_interval_group8{I,J}(:,1),center_cf_interval_group8{I,J}(:,1)]= ...
            hist(cf_interval_group8{I,J},min(cf_interval_group8{I,J}):1:max(cf_interval_group8{I,J}));
        count_cf_interval_group8{I,J}=count_cf_interval_group8{I,J}./sum(count_cf_interval_group8{I,J});
    end
end
for I=1:3
    for J=1:1
        [count_cf_interval_group3{I,J},center_cf_interval_group3{I,J}]=...
            hist(cf_interval_group3{I,J},min(cf_interval_group3{I,J}):1:max(cf_interval_group3{I,J}));
        count_cf_interval_group3{I,J}=count_cf_interval_group3{I,J}./sum(count_cf_interval_group3{I,J});
    end
end
cross_interval.cf_interval_group8=cf_interval_group8;
cross_interval.cf_interval_group3=cf_interval_group3;
cross_interval.count_cf_interval_group8=count_cf_interval_group8;
cross_interval.center_cf_interval_group8=center_cf_interval_group8;
cross_interval.count_cf_interval_group3=count_cf_interval_group3;
cross_interval.center_cf_interval_group3=center_cf_interval_group3;
save('cross_interval.mat','cross_interval');
%%
clear
% load the data
load('D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq\b_MS_SE_STATIC\MS_SE_HY96_RS.mat')
load('dynamic_MS_SE_HY96.mat')

% move between the continuous two window
Move_syn_win = diff(MS_win,1,2);
Move_synE_win = diff(SE_win,1,2);

% calculate the speed for dynamic MS and SE
State_syn_win=MS_win(:,1:100);
Speed_syn_win=abs(Move_syn_win);
State_synE_win=SE_win(:,1:100);
Speed_synE_win=abs(Move_synE_win);
% save
dynamics.move_syn_win=Move_syn_win;
dynamics.move_synE_win=Move_synE_win;
dynamics.state_syn_win=State_syn_win;
dynamics.speed_syn_win=Speed_syn_win;
dynamics.state_synE_win=State_synE_win;
dynamics.speed_synE_win=Speed_synE_win;
save('MS_SE_dynamics_HY96','dynamics')


