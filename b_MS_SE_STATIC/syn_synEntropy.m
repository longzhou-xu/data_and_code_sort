function [MS,SE,KOP] = syn_synEntropy(inputsignal,inputTimeLength,inputCellNumber,numBin)
%this function calculate the synchrony and synchrony entropy of a system
% [syn,syn_entropy,syn_t] = syn_synEntropy(inputsignal,inputTimeLength,inputCellNumber,numBin)
% inputsignal is the system time serise
% inputTimeLength is the number of time points for the system time serise
% inputCellNumber is the number of the system units
% numBin is the hist bin number
% MS is the output return of system's syschronization
% SE is the output return  

    theROITimeCoursesTotal=inputsignal;
    %the hilbert transform in performed to have the analytical signals
    %then have their instaneouse phases
    SIZE=size(theROITimeCoursesTotal);
    Num_vox=SIZE(2);
    if Num_vox~=inputCellNumber
       error('WRONG') 
    end
    Time_Point=SIZE(1);
    if Time_Point~=inputTimeLength
       error('WRONG') 
    end
    for vox=1:Num_vox
        a(:,vox)=hilbert(theROITimeCoursesTotal(:,vox));
        theta(:,vox)=angle(a(:,vox));
    end
    % calculate the Kurumoto order parameteres (r_t) at each time point
    ang=exp(theta*1i);
    for jj=1:1:Time_Point
        KOP(jj)=abs(sum(ang(jj,1:Num_vox)))./Num_vox;
    end    
    MS=mean(KOP);%mean synchrony
    if isnan(MS) == 1
        SE = nan;
    else 
        bins=linspace(min(KOP),max(KOP),numBin);
        b=hist(KOP,bins);
        p=b./Time_Point;
        p(p==0)=[];
        Entropy_rt=0;
        for kk=1:1:length(p)
            Entropy_rt=Entropy_rt-p(kk).*log2(p(kk));
        end
        SE=Entropy_rt;
    end
end

