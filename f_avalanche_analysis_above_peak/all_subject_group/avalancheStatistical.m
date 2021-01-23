function [BinFiringsSerise, AvalancheSizeSerise, DuringTimeSerise, AvaSizeHist, AvaDuraHist] = avalancheStatistical(FiringsSerise, BinSize)
% The ava_size_statistical displays the size of avalanches for the time serise of the firings

% The FiringsSerise is the time serise of the firing or event
% The BinSize the size of the time window

% The AvaSizeHist is the output hist of avalanche size
% The avalanche_size_serise is the size of an avalanches
% The binFirings_serise is the number of event in a time bin
% The duringTime_serise is the duation of an avalanche
%%
% the first step: bin the firing serise; the bin size = inputBinSize
    L = length(FiringsSerise);
    numBin = floor(L / BinSize);
    BinFiringsSerise(:, 1) = sum(reshape(FiringsSerise(1:numBin*BinSize), BinSize, numBin), 1);
    BinFiringsSerise=[0; 0; BinFiringsSerise; 0; 0; ];
% the step two: calculate the duration of avalanche
    event_equal_0 = find(BinFiringsSerise(:,1) == 0);
    DurationN=0;
    for I = 1:length(event_equal_0) - 1
        A = event_equal_0(I+1) - event_equal_0(I);
        if  A > 1
            DurationN = DurationN+1;
            DuringTimeSerise(DurationN,1) = A-1;
        end
    end
% the step three: calculate the size of avalanche
    size_in_a_ava = 0;                      %the size of a avalanche
    J=0;                                  
    AvalancheSizeSerise = [];                          %the serise of avalanche size
    for bin = 1:length(BinFiringsSerise)
        if size_in_a_ava == size_in_a_ava + BinFiringsSerise(bin,1) && size_in_a_ava>0 
           J=J+1;
           AvalancheSizeSerise(J,1) = size_in_a_ava;
           size_in_a_ava = 0;
        elseif size_in_a_ava < size_in_a_ava+BinFiringsSerise(bin,1)
              size_in_a_ava = size_in_a_ava+BinFiringsSerise(bin,1);
        end
    end
% the step four: the distribition of the avalanche size
    AvaSizeHist = zeros(max(AvalancheSizeSerise),1);
    avaNum = length(AvalancheSizeSerise);
    for I = 1:avaNum
        AvaSizeHist(AvalancheSizeSerise(I,1)) = AvaSizeHist(AvalancheSizeSerise(I,1),1)+1;
    end
% the step five: the distribition of the avalanche duration time
    AvaDuraHist = zeros(max(DuringTimeSerise),1);
    for I = 1 : DurationN
        AvaDuraHist(DuringTimeSerise(I,1)) = AvaDuraHist(DuringTimeSerise(I,1),1)+1;
    end
% the step six: test whether correct 
    if avaNum ~= DurationN
        error('wrong')
    end
end

