function [Signals_Z_Normalized]=z_normalize(SignalsOriginal)
% This function is constructed for z-normalizing the original signals.

% The signalsOriginal is the original signals, row is the time, column is the
% ROI.

% The signalsNormalized is the z-normalized signal.

Size = size(SignalsOriginal);
Smean(1,:) = mean(SignalsOriginal, 1);
Smeanmatrix = repmat(Smean,Size(1), 1);
Sstd = std(SignalsOriginal(1:Size(1), :), 1);
Sstdmatrix = repmat(Sstd, Size(1), 1);
Signals_Z_Normalized = (SignalsOriginal(1:Size(1), :)-Smeanmatrix)./Sstdmatrix;

end