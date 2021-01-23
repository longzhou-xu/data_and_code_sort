function [EventSerise, BinEvent, AvalancheSize, DuringTime, AvaSizeHist, AvaDuraHist, BP, BP_other] = avalancheAnalysisBOLD(Signal_normalized,Threshold, InputBinSize)

% calculate the avalanche and its size , duration and branching parameter

% define the event and avalanche
[Event_sum, EventSerise] = eventProduce(Signal_normalized,Threshold);

% calculate the size, duration and their histogram of avalanche,   
[BinEvent, AvalancheSize, DuringTime, AvaSizeHist, AvaDuraHist] = avalancheStatistical(Event_sum, InputBinSize);

% calculate the branchring parameters of each subject
[BP] = branchingParameter(BinEvent);

%%
[BP_other]=branchingParameter_other(BinEvent);

end