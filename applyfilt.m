function [y F] = applyfilt(in,F)
%
% Name: applyfilt
%
% Inputs:
%    in - a scalar, the current filter input
%    F - a struct, representing a digital filter
% Outputs:
%    y - a scalar, the current filter output
%    F - a struct, representing the input filter, updated
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Apply digital filter F to input in, in accordance with 
% the generic filter equation: 
% a(1)*y(i) = b(1)*x(i) + b(2)*x(i-1) + ... + b(n+1)*x(i-n)
%             - a(2)*y(i-1) - ... - a(n+1)*y(i-n)
%

% Determine output on basis of recent inputs (including 'in')
% recent output and the filter coefficients
y = (F.b'*[in; F.in] - F.a(2:end)'*F.out)/F.a(1);

% Update recent inputs
%    'in' is now the most recent
%    least recent is discarded
F.in = [in; F.in(1:end-1)];

% Update recent outputs
%    'y' is now the most recent
%    least recent is discarded
F.out = [y; F.out(1:end-1)];

return
%eof