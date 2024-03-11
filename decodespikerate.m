function x_t_hat = decodespikerate(At,PHI)
%
% Name: decodespikerate
%
% Inputs:
%    At - p-by-n vector of LIF temporal activation functions,
%         columns corresponding to the n individual neurons in the pop'ln,
%         rows corresponding to the p samples in the temporal activation function.
%    PHI - n-by-1 vector of decoders,
%          rows corresponding to the n individual neurons in the pop'ln
% Outputs:
%    x_t_hat - p-by-1 estimate of the time-varying value encoded by the 
%               population of neurons      
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Estimate the value encoded by a population of neurons on
%              the basis of their spiking rates recorded at a common
%              physical value, as well as their determined decoders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are at least three viable solutions to this function,
% depending on which type of Matlab syntax is employed:

p = size(At,1);
n = size(At,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #1: Using Iteration Over Neurons

% %%% Determine Estimate of Physical Value x
% x_t_hat = zeros(p,1);
% for itor = 1:n
%     x_t_hat = x_t_hat + At(:,itor)*PHI(itor);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #2: Using Iteration Over Samples

%%% Determine Estimate of Physical Value x
% x_t_hat = zeros(p,1);
% for itor = 1:p
%     x_t_hat = x_t_hat + sum(At(itor,:).*PHI');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #3: Using Vector Operations - e.g., the inner product

% %%% Determine Estimate of Physical Value x
x_t_hat = At*PHI;

return
%eof