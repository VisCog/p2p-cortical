function [response, blanks] = ecogSubtractBlanks(response, stimulus)

% subtract out baseline: which stimulus frames have no contrast?
blanks = mean(stimulus, 2) < 10e-6;

if sum(blanks) == 0
    % do nothing - there are no blanks
    
else
    % subtract the mean of each run from that run's time series
    response = response - mean(response(blanks));
end