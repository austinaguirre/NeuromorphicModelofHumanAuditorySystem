function N = updatelifn(x,N,Fs)
%
% Name: updatelifn
%
% Inputs:
%    x - a (scalar) physical value, to be encoded
%    N - a struct, representing an LIF neuron to be updated
%    Fs - a (scalar) sampling rate, in Hz
% Outputs:
%    N - a struct, representing the updated LIF neuron
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Update the membrane voltage an LIF neuron, 
%               on the basis of its parameters and a physical value 

tau = 1/Fs; % simulation time step in seconds

if N.state == 1 % Super-Threshold Neuron Behavior
    
    if N.refc==0 % If refractory period "clock" has been reset (==0)
        N.V = 1; % Set membrane voltage to '1' (spike)
    else         % Else
        N.V = 0; % Set membrane voltage to '0' (refractory period voltage)
    end
    
    if N.refc < N.tau_ref     % If refractory period "clock" is less than absolute refractory period
        N.refc = N.refc + tau;% Add time to the "clock" equal to one simulation time step
    else                      % Else (i.e., the refractory period is over)
        N.state = -1;         % Change state to '-1' (sub-threshold behavior)
        N.refc = 0;           % Reset the refractory period "clock" (=0)
    end
    
else % Sub-Threshold Neuron Behavior (State == '-1')
    
    J = N.alpha*N.phi*x + N.J_bias; % Determine soma current (see 'ratelifn.m')
    
    Vd = -(1/N.tau_RC)*(N.V-J*N.R); % Determine change in membrane voltage (see Eliasmith 4.4)
    
    N.V = N.V + Vd*tau; % Update membrane voltage using Euler integration
    
    if N.V >= N.Vth % If membrane voltage is above threshold
        N.state = 1; % Change state to '1' (super-threshold behavior)
    end
    
end
            

return
%eof