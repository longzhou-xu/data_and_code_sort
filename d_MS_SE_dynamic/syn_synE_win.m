function [Mean_synchrony_win,Synchrony_entropy_win] = syn_synE_win(KOP,T,Bin,Slide_size,Num_entropy)
%[syn_win,synE_win] = syn_synE_win(KOP,T,bin,slide)
%calculate the mean synchrony and synchrony entropy in slide window
Size=size(KOP);
if T ~= Size(1)
    error('the number of timepoint wrong');
end
Win_end = Bin:Slide_size:Size(1); % the end index of every window
Win_begin = Win_end-Bin+1;% the begin index of every window

Mean_synchrony_win = zeros(Size(2),length(Win_end));
Synchrony_entropy_win = zeros(Size(2),length(Win_end));

for Win_num = 1:length(Win_end) % traverse all the window
    
    KOP_win(1:Bin,1:Size(2)) = KOP(Win_begin(Win_num):Win_end(Win_num),:); % the KOP in the Win_num~th window 
    
    Mean_synchrony_win(:,Win_num) = mean(KOP_win,1); % the mean synchrony in the Win_num~th window
    
    for Subject = 1:Size(2) % traverse all the subject to calculate the synchrony entropy
        if isnan(Mean_synchrony_win(Subject,Win_num)) == 1
           Synchrony_entropy_win(Subject,Win_num) = nan;
        else
           Bins=linspace(min(KOP_win(:,Subject)),max(KOP_win(:,Subject)),Num_entropy);
           B=hist(KOP_win(:,Subject),Bins);
           P=B./sum(B);
           P(P==0)=[];
           Entropy_rt=0;
           for kk=1:1:length(P)
               Entropy_rt=Entropy_rt-P(kk).*log2(P(kk));
           end
           Synchrony_entropy_win(Subject,Win_num)=Entropy_rt;
        end
        
        clear P Entropy_rt
    end
    clear KOP_win
end
end

