function [fitresult, fitX, fitY, gof, output] = f_fit_poly1(X, Y)
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


%% main
    [xData, yData] = prepareCurveData( X, Y);

    % Set up fittype and options.
    ft = fittype( 'poly1' );

    % Fit model to data.
    [fitresult, gof, output] = fit( xData, yData, ft );

    p1 = fitresult.p1; 
    p2 = fitresult.p2; 

    d = (max(X) - min(X))/100;

    fitX = min(X)-2*d : d : max(X)+2*d;

    fitY = p1 .* fitX + p2;

end


