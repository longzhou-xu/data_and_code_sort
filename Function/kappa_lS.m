function [Kappa] = kappa_lS(ava_size,X_min,X_max,bin_num)
%calculate the kappa

%% the step one 
% from the sample to calculate the cumulative density function: F
% the size range
if X_max > max(ava_size)
    X_max = max(ava_size);
end
betak = logspace(log10(X_min),log10(X_max),bin_num);
ava_size(ava_size>X_max) = [];
ava_size(ava_size<X_min) = [];
[size_hist] = hist(ava_size,betak);
mava = size_hist ./ sum(size_hist);
  
%input is the avalanches size serise
    L = X_max;
    size_hist = hist(ava_size,X_min:X_max);
    size_hist_p = size_hist./sum(size_hist);
    syms x;
    equ = (sqrt(X_max) / (sqrt(X_max)-sqrt(x))) * (1 - sqrt(x) ./ sqrt(X_min)) == size_hist_p(1);
    m_FNA = solve(equ,x);
    l = eval(m_FNA);
    Z = sqrt(L)./(sqrt(L)-sqrt(l));
    
    kappa = 0;

    for j = 1:length(betak)
        a = sqrt(betak(j));
        b = sqrt(l);
        FNA = Z * (1 - b./a);
        F = sum(mava(1:j));
        kappa = kappa + FNA - F;
        FNA_OUT(j)=FNA;
        F_OUT(j)=F;
    end
    
    kappa = 1 + kappa / bin_num;
    
    Kappa.kappa = kappa;
    Kappa.max_size = X_max;
    Kappa.min_size = X_min;
    Kappa.L = L;
    Kappa.betak = betak;
    Kappa.F_data = F_OUT;
    Kappa.F_na = FNA_OUT;
    
    figure('Name','kappa','Color','w');
    F1 = plot(1:j, FNA_out,...
        'LineStyle', '--', ...
        'Color', 'b');% plot the theoritical CDF
    hold on
    F2 = plot(1:j, F_OUT, ...
        'LineStyle', '-',...
        'Color', 'r');% plot the data CDF 
    lenged([F1,F2], {'FNA','F'}, ...
        'Interpreter', 'latex', ...
        'FontName', 'Arial', ...
        'FontSize', 20, ...
        'EdgeColor', [1,1,1]);
    
end