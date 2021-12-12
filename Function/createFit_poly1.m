function [fitresult, gof, output] = createFit_poly1(syn_rest_include, corrFCmean_include)
%CREATEFIT(SYN_REST_INCLUDE,CORRFCMEAN_INCLUDE)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : syn_rest_include
%      Y Output: corrFCmean_include
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 15-Dec-2019 22:26:03 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( syn_rest_include, corrFCmean_include );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'corrFCmean_include vs. syn_rest_include', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'syn_rest_include', 'Interpreter', 'none' );
% ylabel( 'corrFCmean_include', 'Interpreter', 'none' );
% grid on
end


