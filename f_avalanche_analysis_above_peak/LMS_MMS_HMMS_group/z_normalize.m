function [signalsNormalized]=z_normalize(signalsOriginal)
% This function is constructed for z-normalizing the original signals.

% The signalsOriginal is the original signals, row is the time, column is the
% ROI.

% The signalsNormalized is the z-normalized signal.

Size = size(signalsOriginal);
Smean(1,:) = mean(signalsOriginal, 1);
Smeanmatrix = repmat(Smean,Size(1), 1);
Sstd = std(signalsOriginal(1:Size(1), :), 1);
Sstdmatrix = repmat(Sstd, Size(1), 1);
signalsNormalized = (signalsOriginal(1:Size(1), :)-Smeanmatrix)./Sstdmatrix;

end