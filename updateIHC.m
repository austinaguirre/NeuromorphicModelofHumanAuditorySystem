function IHC = updateIHC(input_signal, IHC)
% Update the IHC model based on an input signal
% Inputs:
%   input_signal - A scalar representing the input signal to the IHC
%   IHC - A struct representing the IHC model
% Outputs:
%   IHC - The updated IHC model
%     'V' - IHC membrane voltage

tau = 1 / IHC.fs; % Simulation time step in seconds

J = IHC.alpha * IHC.phi * input_signal + IHC.J_bias; % Calculate the current

% Update the membrane voltage of the IHC using Euler integration
Vd = -(1 / IHC.tau_RC) * (IHC.V - J * IHC.R);
IHC.V = IHC.V + Vd * tau;

return
%eof