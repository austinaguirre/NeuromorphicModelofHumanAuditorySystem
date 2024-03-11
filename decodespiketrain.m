function x_hat_j = decodespiketrain(N,D,Fs)
%
% Name: decodespiketrain
%
% Inputs:
%    N - A struct, representing a LIF neuron
%    D - a n-by-1 vector/signal, representing a spike train recording
%    Fs - a (scalar) sampling rate, in Hz
% Outputs:
%    x_hat_j - a n-by-1 vector/signal, representing the temporal activation
%           function of the input neuron. 
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Calculate a partial estimate, attributed to one neuron,
%              of the value encoded by a population of neurons, on
%              the basis of their spiking rates recorded at a common
%              physical value, as well as their determined decoders

t = linspace(0,length(D)/Fs,length(D))'; % Create vector of time stamps (sec)

h_t = exp(-t/(1*N.tau_RC)); % Create post-synaptic current (PSC) vector

x_hat_j = conv(h_t,D); % Convolve decoder-weighted PSC with spike train

return
%eof