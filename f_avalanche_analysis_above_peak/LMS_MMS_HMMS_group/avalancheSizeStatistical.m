function [AvaSizeHist, avalanche_size_serise, binFirings_serise, duringTime_serise] = avalancheSizeStatistical(inputFiringsSerise, BinSize)
% The ava_size_statistical displays the size of avalanches for the time serise of the firings

% The FiringsSerise is the time serise of the firing or event
% The BinSize the size of the time window

% The AvaSizeHist is the output hist of avalanche sizw
% The avalanche_size_serise is the size of an avalanches
% The binFirings_serise is the number of event in a time bin
% The duringTime_serise is the duation of an

% the first step: bin the firing serise; the bin size = inputBinSize
    l = length(inputFiringsSerise);
    numBin = floor(l / BinSize);
    binFirings_serise(:, 1) = sum(reshape(inputFiringsSerise(1:BinSize*numBin), numBin, BinSize), 2);
    binFirings_serise=[0; 0; binFirings_serise; 0; 0; ];
% the step two: calculate the duration of avalanche
    event_equal_0 = find(binFirings_serise(:,1) == 0);
    K=0;
    for I = 1:length(event_equal_0) - 1
        A = event_equal_0(I+1) - event_equal_0(I);
        if  A > 1
            K = K+1;
            duringTime_serise(K,1) = A-1;
        end
    end
% the step three: calculate the size of avalanche
    size_in_a_ava = 0;                      %the size of a avalanche
    J=0;                                  
    avalanche_size_serise = [];                          %the serise of avalanche size
    for bin = 1:length(binFirings_serise)
        if size_in_a_ava == size_in_a_ava + binFirings_serise(bin,1) && size_in_a_ava>0 
           J=J+1;
           avalanche_size_serise(J,1) = size_in_a_ava;
           size_in_a_ava = 0;
        elseif size_in_a_ava < size_in_a_ava+binFirings_serise(bin,1)
              size_in_a_ava = size_in_a_ava+binFirings_serise(bin,1);
        end
    end
% the step four: the distribition of the avalanche size
    AvaSizeHist = zeros(max(avalanche_size_serise),1);
    avaNum = length(avalanche_size_serise);
    for I = 1:avaNum
        AvaSizeHist(avalanche_size_serise(I,1)) = AvaSizeHist(avalanche_size_serise(I,1),1)+1;
    end
% the step five: test whether correct 
    if avaNum ~= K
        error('wrong')
    end
end

