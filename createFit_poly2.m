function [fitresult, gof, output] = createFit_poly2(syn_rest_include, connect_number_entropy1_include)
%CREATEFIT(SYN_REST_INCLUDE,CONNECT_NUMBER_ENTROPY1_INCLUDE)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : syn_rest_include
%      Y Output: connect_number_entropy1_include
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( syn_rest_include, connect_number_entropy1_include );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'connect_number_entropy1_include vs. syn_rest_include', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'syn_rest_include', 'Interpreter', 'none' );
% ylabel( 'connect_number_entropy1_include', 'Interpreter', 'none' );
% grid on
% 
end
