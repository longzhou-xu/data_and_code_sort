function [mean_syn, std_syn, cv_syn, entropy_syn,sample_failed] =  kop2msse(kuramoto_order_parameter,bins_num)
%% notes:
%% input:
%% output:
%% main
        Time_Point = length(kuramoto_order_parameter);
        % mean of synchrony
        mean_syn = mean(kuramoto_order_parameter);
        % std of synchrony 
        std_syn = std(kuramoto_order_parameter);
        % cv of synchrony
        cv_syn = std_syn/mean_syn;
        
        % [0, 1] normalize
        kop_max = max(max(kuramoto_order_parameter));
        kop_min = min(min(kuramoto_order_parameter));
        kop_nor = (kuramoto_order_parameter - kop_min)./(kop_max - kop_min); 
        
        % calculate the synchrony entropy
        bins = linspace(0, 1, bins_num);
        bb = hist(kop_nor, bins);
        sample_failed = length(find(bb == 0));
        
        prob = bb./Time_Point;
        prob(prob == 0) = [];
        
        entropy_syn = 0;
        for kk = 1:1:length(prob)
            entropy_syn = entropy_syn - prob(kk) .* log2(prob(kk));
        end
end