function a = ratelifn(N,x)
%
% Name: ratelifn
%
% Inputs:
%    N - A struct, representing a LIF neuron
%    x - A scalar value, to be encoded 
% Outputs:
%    a - LIF neuron spiking rate (spikes/sec)
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Calculate the spiking rate of an LIF neuron
%               given a specified value

% Determine soma currents
J_i = N.alpha*N.phi*x + N.J_bias;

% Determine spiking rates
if J_i > N.J_th
    G_i = 1./(N.tau_ref - N.tau_RC*log(1-(N.J_th./J_i)));
else
    G_i = 0;
end

a = G_i;

return
%eof