function [ava] = detect_avalanches(raster, time_length, timebinSize)
%% notes
%
%% input
%
%% output
%
%% main function
% sum the event
    event_sum = sum(raster,1);
% bin the event serise
    t_len = length(event_sum);
    if time_length ~= t_len
        error('time dimension wrong')
    end
    num_bins = floor(t_len / timebinSize);
    event_bin(1,:) = sum(reshape(event_sum(1 : num_bins*timebinSize), num_bins, timebinSize), 2);
    event_bin = [0, event_bin, 0];
% detect the avalanche and calculate the property of each avalanche
    iZeroBin = find(event_bin == 0);
    iav = 0; % the index of avalanche  
    for tt = 1 : length(iZeroBin) - 1
        AA = iZeroBin(tt+1) - iZeroBin(tt);
        if  AA > 1
            iav = iav + 1; 
            
            ava.begin_time(iav) = iZeroBin(tt) + 1 - 1;% the begin time point of each avalanche
            ava.end_time(iav) = iZeroBin(tt + 1) - 1 - 1;% the end time point of each avalanche
            
            ava.shape{iav} = event_bin(iZeroBin(tt) + 1 : iZeroBin(tt+1) - 1);
            
            % calculate the duration of each avalanche
            ava.duration(iav) = ava.end_time(iav) - ava.begin_time(iav) + 1;
            
            % calculate the size of each avalanche
            ava.size(iav) = sum(ava.shape{iav});
            
            % convert the raster to fingerprint
            clear raster_iav 
            raster_iav = raster(:, iZeroBin(tt) + 1 - 1 : iZeroBin(tt+1) - 1 - 1);
            [row, col] = find(raster_iav == 1);
            ava.fingerprint(iav).time = col + ava.begin_time(iav)  - 1;
            ava.fingerprint(iav).node = row;
            
            % check the possible wrong
            if sum(ava.shape{iav} == 0) > 1 
                error('wrong algorithm');
            end
            
        end
    end
end

