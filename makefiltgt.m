function F = makefiltgt(fc,bw,Fs)
%
% Name: makefiltgt
%
% Inputs:
%    fc - Center frequency (Hz)
%    bw - Filter bandwidth (Hz)
%    Fs - Sampling frequency (Hz)
% Outputs:
%    F - a struct, representing a 2nd-order gammatone digital filter
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Create a 2nd-order IIR digital filter
%              and populate its coefficients in accordance with a
%              gammatone filter with the specified cutoff frequency
%              and bandwidth at sampling rate Fs
%

% Make a 2nd-order digital filter
F = makefilt(2);

% Determine filter cutoff frequency and bandwidth in radians/sample
theta = (pi*fc)/(Fs/2);
beta = (pi*bw)/(Fs/2);
r = exp(-beta);

% Determine additional parameters
j = sqrt(-1);
alpha = r*cos(theta);

% Calculate feedback coefficients
% Store them in the filter struct
F.a(1) = 1;
F.a(2) = -2*alpha;
F.a(3) = r^2;

% Calculate feedforward coefficients
% Store them in the filter struct
F.b(1) = abs(...
    (1+F.a(2)*exp(-j*theta) + F.a(3)*exp(-j*2*theta))/...
    (1-r*cos(theta)*exp(j*theta))...
    );

F.b(2) = -alpha*F.b(1);
F.b(3) = 0;

return
%eof