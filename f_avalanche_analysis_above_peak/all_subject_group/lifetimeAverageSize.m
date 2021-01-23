function [Lifetime, AverageSize]=lifetimeAverageSize(AvaDuration, AvaSize)
% The lifetimeAverageSize is of calculating the average size for same 
% AvaDuration: duration of avalanche.
% AvaSize: 
Max_lifetime = max(AvaDuration);
Size_sort = cell(Max_lifetime, 1);
for L_AvaD = 1:length(AvaDuration)
    N = Size_sort{AvaDuration(L_AvaD, 1), 1};
    N = [N; AvaSize(L_AvaD, 1)];
    Size_sort{AvaDuration(L_AvaD,1), 1}=N;
end
AverageSize = zeros(length(Size_sort), 1);
for M=1:length(Size_sort)
    AverageSize(M,1)=mean(Size_sort{M, 1});
end
Lifetime = find(AverageSize>0);
AverageSize(isnan(AverageSize)) = [];
end