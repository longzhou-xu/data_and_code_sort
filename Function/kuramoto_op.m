function [kop_r, kop_complex] = kuramoto_op(SIGNALS)
%% document
% [kura_order_pm, kop_complex] = kuramoto_op(signals, time_len, node_num)
% This function calculate the Kuramoto order parameter from the signals.
% First, we use the Hilbert transform to obtain the analytic signals (complex number) for data.
% Second, we calculate the instaneouse phase from the analytic signals.
% Finaly, we calculate the Kuramoto order parameter using the transient phase.

% Input
% SIGANLS are the siganls data from one system. 
% the row of SIGNALS is index of elements.
% the column of SIGNALS is index of time. 

% output

%% main function
    [e_num, t_len] = size(SIGNALS);
    
    % obtain the instaneouse phase
    for NN = 1 : e_num
        
        a_s(NN,:)=hilbert(SIGNALS(NN,:)); %analytical signals

        theta(NN,:)=angle(a_s(NN,:));% instaneouse phases
    end
    
    phase = exp(theta * 1i); % the unit vector of phase
    
    % calculate the Kurumoto order parameteres at each time point
    kop_complex = sum(phase, 1);
    
    kop_r = abs(kop_complex) / e_num;
end

