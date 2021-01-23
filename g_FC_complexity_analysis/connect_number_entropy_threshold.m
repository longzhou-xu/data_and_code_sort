function [CNE_subject,CNE_ROI] = connect_number_entropy_threshold(Bold,ROI_number,threshold,bin_number,hist_number)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
SIZE=size(Bold);
if SIZE(2)~=ROI_number
    error('the input is wrong')
end
%% 计算滑动窗口的FC矩阵
win_head = floor(linspace(1,SIZE(1),bin_number));  % 把整个时间序列近似线性等分
corrMatrix = zeros(length(win_head)-1, ROI_number, ROI_number);
for i = 1 : length(win_head) - 1   %对时间窗口的个数循环
    win_b = win_head(i);  %时间窗口的起始点位置
    win_e = win_head(i + 1); %时间窗口的结束点位置
    corrMatrix(i, 1 : ROI_number, 1 : ROI_number) = corr(Bold(win_b : win_e, 1 : ROI_number));%此3维矩阵为当前被试所有时间窗口的关联矩阵
    corrMatrix(corrMatrix >= 0.9999) = 0;    %把关联矩阵对角线上的1换成0
end
%% 对FC矩阵进行二值化建立网络
    corrMatrix = abs(corrMatrix);  %对关联矩阵取绝对值
    corrMatrix(corrMatrix >= threshold) = 1; corrMatrix(corrMatrix <= threshold) = 0;%对关联矩阵中大于阈值的取为1，小于的取为0
%% 对每一个脑区计算CNE 
    for i = 1 : ROI_number    %对脑区循环
        n1(1 : length(corrMatrix(:, 1, 1)), 1 : ROI_number) = corrMatrix(: , i, 1 : ROI_number);
        n2 = sum(n1,2);  %每个时间窗口中，当前脑区和其他脑区关联的个数
        p = hist(n2,linspace(1,ROI_number,hist_number))./length(n2);  %除以总数求得概率
        p(p == 0) = [];   %去掉概率为0的
        p = -p.*log2(p); %此行及下行 为熵的公式
        CNE_ROI(i,1) = sum(p);
        clear p
%         fprintf('%d: %d\n',i,sum(n2)); %计算过程中输出每个脑区的连接个数，以判断阈值过大或过小
    end 
 CNE_subject = mean(CNE_ROI);   %被试的平均CNE
end

