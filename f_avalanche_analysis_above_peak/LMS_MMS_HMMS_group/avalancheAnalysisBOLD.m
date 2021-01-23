function [avalancheSize, binFiring, duringTime, branch] = avalancheAnalysisBOLD(signal_normalized,threshold, inputBinSize)

% calculate the avalanche and its size , duration and branching parameter

% define the event and avalanche
[event_sum, event] = BOLDeventProduce(signal_normalized,threshold);

% calculate the size, duration and their histogram of avalanche,   
[SizeH, avalancheSize, binFiring, duringTime] = avalancheSizeStatistical(event_sum, inputBinSize);

% calculate the branchring parameters of each subject
[branch] = branchingParameter(binFiring);

end