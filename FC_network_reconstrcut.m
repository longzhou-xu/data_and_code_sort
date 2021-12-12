clc;clear;close all
%
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
%
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\b_MS_SE_STATIC\MS_SE_Z1024_RS.mat'], 'subject')
mkdir('FC\FC_binary')
load(['D:\criticality_cognitive\project_synchrony_avalanche_fluid_Iq',...
    '\g_FC_complexity_analysis\static_FC_complexity.mat'], 'FC')
for S=1:295
    disp(S)
    clear Cij
    Cij=abs(FC(:,:,S));%Cij为原始功能连接矩阵，Aij为二值化之后的功能连接矩阵
    Cij=Cij-diag(diag(Cij));%对角元素归零（不考虑自身相关）
    %通过控制连接密度来确定二值化阈值
    %连接密度0.01：0.01：0.99
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
    save(['FC\FC_binary\sub_',num2str(S)],'FC_binary'); 
    save('FC\subject.mat','subject');
    clear FC_binary
end