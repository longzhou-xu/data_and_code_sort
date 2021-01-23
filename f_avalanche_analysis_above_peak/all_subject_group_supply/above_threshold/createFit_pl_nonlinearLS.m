function [fitresult, gof, output] = createFit_pl_nonlinearLS(x, y)
%CREATEFIT(SIZE_RANGE,AVA_SIZE_FIT_P)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : size_range
%      Y Output: ava_size_fit_p
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.
%  由 MATLAB 于 05-Dec-2019 10:16:05 自动生成
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Robust = 'LAR';
% opts.Algorithm = 'Levenberg-Marquardt';
opts.StartPoint = [0.43366588132904 -1.55666393334776];


% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );

end


