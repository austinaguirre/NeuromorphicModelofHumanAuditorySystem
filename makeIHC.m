function IHC = makeIHC(alpha, phi, J_bias, fs)
% Create an Inner Hair Cell (IHC) model
% Inputs:
%   alpha - Scaling parameter for IHC
%   phi - Encoding weight for IHC
%   J_bias - Bias or background current for IHC
%   fs - Sampling rate in Hz
% Outputs:
%   IHC - A struct representing the IHC
%     'V' - IHC membrane voltage
%     'alpha' - Scaling parameter
%     'phi' - Encoding weight
%     'J_bias' - Bias or background current
%     'fs' - Sampling rate
%     'tau_RC' - RC time constant
%     'R' - Leak resistance

IHC = struct(...
    'V', 0, ...
    'alpha', alpha, ...
    'phi', phi, ...
    'J_bias', J_bias, ...
    'fs', fs, ...
    'R',10,...
    'tau_RC',0.0001...
    );

return
%eof