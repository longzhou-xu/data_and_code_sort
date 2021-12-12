function [H, zeroN] = entropy(data,bins)
%% note
%% input
%% output
%% main 
    b = hist(data, bins);
    p = b ./ length(data);
    zeroN = length(find(p == 0));
    p(p == 0) = [];
    H = 0;
    for binN = 1 : 1 : length(p)
        H = H - p(binN) .* log2(p(binN));
    end
end