function N = makelifn(alpha,phi,J_bias)

N = struct(...
    'state',-1,... % spike state: -1,subthreshold; 1,superthreshold
    'V',0,... % LIF membrane voltage
    'refc',0,... % refractory period "clock"
    'alpha',alpha,... % scaling parameter
    'phi',phi,... % encoding weight
    'J_bias',J_bias,... % "bias" or "background" current
    'J_th',0,... % current threshold
    'tau_RC',0.0001,... % RC time constant
    'tau_ref',0.00001,... % absolute refractory period
    'R',10,... % leak resistance
    'Vth',0.62... % LIF voltage threshold
    );

N.J_th = N.Vth/N.R;

return
%eof