function [EventSumSerise,BOLDeventSerise]=eventProduce(BOLDsignals,Threshold)
% The BOLDeventProduce define the event of avalanche of BOLD siganls.

% The BOLDsignals is the z-normalized signals.
% The threshold is the BOLD threshold for defining the envent.

% The eventSumSerise is the sum of event across ROI.
% The BOLDeventSerise records the event and its position and time.
% The events defined as the above threshold. 
%%
    Size = size(BOLDsignals);
    BOLDeventSerise = zeros(Size(1), Size(2));
    for roi=1 : Size(2)
        for T=1 : Size(1)-2
            if  BOLDsignals(T, roi) >= Threshold && BOLDsignals(T+1, roi)>= Threshold ...
                    && BOLDsignals(T+2, roi) >= Threshold
                if BOLDsignals(T, roi) <= BOLDsignals(T+1, roi) && ...
                                 BOLDsignals(T+1, roi) >= BOLDsignals(T+2, roi)
                   BOLDeventSerise(T+1, roi) = 1;
                end
            end
        end
    end
    EventSumSerise(:,1) = sum(BOLDeventSerise,2);  
end