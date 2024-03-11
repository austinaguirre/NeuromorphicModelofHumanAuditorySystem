function [A X] = characterizelifn(N,minval,maxval,res,Fs)
%
% Name: characterizelifn
%
% Inputs:
%    N - A struct, representing a LIF neuron
%    minval - The minimum relevant value, to be encoded
%    maxval - The maximum relevant value, to be encoded
%    res - The number of relevant values (evenly-spaced) 
%           to be encoded over the range specified by minval and maxval
%    Fs - a (scalar) sampling rate, in Hz
% Outputs:
%    A - res-by-1 vector of LIF neuron spiking rates (spikes/sec),
%         corresponding to evenly-spaced values between minval and maxval
%    X - a res-by-1 vector of values, to be encoded
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Estimate the spiking rate of an LIF neuron 
%               by directly observing spiking behavior, given
%               a specified range of values and
%               a number of evenly-spaced values over that range

nsamp = floor(Fs./10); % Calculate the number of samples required 
                       % for good estimate of spiking rate

X = linspace(minval,maxval,res)'; % Create a vector of relevant values

A = zeros(res,1); % Initialize vector of probed spiking rates

% Iterate over relevant values
for itor = 1:res
    
    counter = 0; % Initialize spike counter (=0)
    
    % Iterate over the required number of samples
    for jtor = 1:nsamp
        
        N = updatelifn(X(itor),N,Fs); % Update the neuron
        
        if N.V == 1                % If membrane voltage is '1' (spike)
            counter = counter + 1; % Increment the spike counter
        end
        
    end
    
    A(itor) = counter*(Fs/nsamp); % Record spike total in spikes/second
    
end

return
%eof