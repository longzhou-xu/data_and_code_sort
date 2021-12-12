clc;clear;close all

den=(1:(95*96/2))./(95*96/2);
I=0;
Den_index=zeros(80,2);
for Den=0.01:0.01:0.99
    I=I+1;
    a=find(den>Den);
    Den_index(I,1)=Den;
    if a(1)/(95*96/2)-Den < (a(1)-1)/(95*96/2)-Den
        Den_index(I,2)=a(1);

    else
        Den_index(I,2)=a(1)-1;
    end
end
load('FC_complexity_flexibility\static_FC_complexity_original_signal.mat', 'FC')
filedir=dir('FC_complexity_flexibility\FC_win_original_signal');
filedir=filedir(3:end,1);
mkdir('FC_SC_similarity\binary_FC_original_signal');
for Subject=1:length(filedir)
    disp(Subject)
    clear Cij
    Cij=abs(FC(:,:,Subject));
    Cij=Cij-diag(diag(Cij));
    %link density : 0.01£º0.01£º0.99
    I=0;
    FC_seise=zeros(95*96/2,1);
    for roi1=1:96
        for roi2=roi1+1:96
            I=I+1;
            FC_seise(I,1)=Cij(roi1,roi2);
        end
    end
    [FC_serise_sort,FC_serise_sort_index]=sort(FC_seise,'descend');
    J=0;
    FC_binary=zeros(96,96,99);
    for Den = 0.01:0.01:0.99
        J=J+1;
        Aij=Cij;
        Aij(Aij<FC_serise_sort(Den_index(J,2)))=0;
        FC_binary(:,:,J)=Aij;
    end
    save(['FC_SC_similarity\binary_FC_original_signal\',filedir(Subject,1).name(1:6)],'FC_binary');   
end

