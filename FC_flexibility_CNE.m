%% FC felxibility CNE
clc; clear; close all
CNE_subject=zeros(295,1,1);
CNE_ROI=zeros(295,96,1,1);
for subject=1:295
    subject
    load(['HY_96\sub', num2str(subject), '.mat']);

    Cmean(1,:) = mean(rest_HY96_ROI(1:1200,:),1);
    Cmeanmatrix = repmat(Cmean,1200,1);
    Cstd = std(rest_HY96_ROI(1:1200,:),1);
    Cstdmatrix = repmat(Cstd,1200,1);
    rest_C = (rest_HY96_ROI(1:1200,:)-Cmeanmatrix)./Cstdmatrix;
    
    I=0;
    for threshold = 0.3 : 0.02 : 0.3 % 0.12 : 0.02 : 0.48
        I=I+1;
        J=0;
        for win_number = 60 : 2 : 60 % 40 : 2 : 60
            J=J+1;
            BOLD = rest_C; ROI_number=96; hist_number=30;
            [CNE_subject(subject, I, J),CNE_ROI(subject, :, I, J)] = ...
                connect_number_entropy_threshold(BOLD,ROI_number,...
                threshold,win_number,hist_number);
        end
    end
    
end
save('flexibility', 'CNE_subject', 'CNE_ROI');
