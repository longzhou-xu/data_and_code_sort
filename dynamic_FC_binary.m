clc; clear; close all
% the density 
den=(1:(95*96/2))./(95*96/2);
I=0;
Den_index=zeros(80,2);
for Den=0.01:0.01:0.8
    I=I+1;
    a=find(den>Den);
    Den_index(I,1)=Den;
    if a(1)/(95*96/2)-Den < (a(1)-1)/(95*96/2)-Den
        Den_index(I,2)=a(1);

    else
        Den_index(I,2)=a(1)-1;
    end
end

for S=1:295
    load(['FC_win\sub',num2str(S),'.mat'])
    disp(S)
    for  win=1:101
         clear Cij
         Cij=abs(FC_win(:,:,win));
         Cij=Cij-diag(diag(Cij));

         %link density : 0.01£º0.01£º0.8
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
         FC_binary=zeros(96,96,80);
         for Den = 0.01:0.01:0.8
             J=J+1;
             Aij=Cij;
             Aij(Aij<FC_serise_sort(Den_index(J,2)))=0;
             FC_binary(:,:,J)=Aij;
         end
         mkdir(['FC_win_binary\sub',num2str(S),'\win',num2str(win)])
         save(['FC_win_binary\sub',num2str(S),'\win',num2str(win),'\binary_FC'],'FC_binary');   
    end
end
