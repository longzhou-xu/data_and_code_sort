function [raster,fingerprint]=detect_events(signals,threshold, time_length, node_num)
%% note

%% input

%% output

%% main function
    t_len = length(signals(1,:));
    n_num = length(signals(:,1));
    % check the dimension
    if time_length ~= t_len
        error('time dimension wrong')
    end
    if node_num ~= n_num
        error('node dimension wrong')
    end
    % check the if the signals are zscore normalized
    if sum(std(signals')) ~= n_num || sum(mean(signals, 2)) > 0.0001
        error('signals without normalizing')
    end
    % convert the signals to raster
    raster = zeros(n_num, t_len);
    for NN = 1 : n_num
        for TT = 1 : t_len-2
            if signals(NN,TT) >= threshold && signals(NN,TT+1) >= threshold ...
                    && signals(NN,TT+2) >= threshold
                
                if  signals(NN, TT) <= signals(NN,TT+1) ...
                        && signals(NN, TT+1) >= signals(NN,TT+2)
                    raster(NN, TT+1) = 1;
                end
                
            end
        end
    end
    % convert the raster to fingerprint
    [row,col] = find(raster == 1);
    fingerprint.time = col;
    fingerprint.node = row;
end