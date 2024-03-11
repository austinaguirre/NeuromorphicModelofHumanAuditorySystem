function PHI = determinedecoders(A,X)
%
% Name: determinedecoders
%
% Inputs:
%    A - p-by-n matrix of LIF neuron spiking rates (spikes/sec),
%         rows corresponding to the p physical values contained in X,
%         columns corresponding to the n individual neurons in the pop'ln.
%    X - p-by-1 vector of physical values, 
%         encoded in spiking rates in A (one for each row)
% Outputs:
%    PHI - n-by-1 vector of decoders
%          rows corresponding to the p individual neurons in the pop'ln
%
% Created by: Adam C. Lammert (2022)
% Author: ??? (you)
%
% Description: Determine decoders for a population of neurons, on the
%              basis of their activation functions - i.e., their 
%              spiking rates measured at a range of physical values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are at least three viable solutions to determine Gamma and Upsilon,
% depending on which type of Matlab syntax is employed:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #1: Using Iteration - e.g., for-loops

% p = size(A,1); % number of physical values
% n = size(A,2); % number of neurons/decoders
% 
% %%% Determine Upsilon
% Upsilon = zeros(n,1);
% for ntor = 1:n
%     for ptor = 1:p
%         Upsilon(ntor) = Upsilon(ntor) + X(ptor)*A(ptor,ntor);
%     end
% end
% 
% %%% Determine Gamma
% Gamma = zeros(n,n);
% for itor = 1:n
%     for jtor = 1:n
%         for ptor = 1:p
%             Gamma(itor,jtor) = Gamma(itor,jtor) + A(ptor,itor)*A(ptor,jtor);
%         end
%     end
% end
% 
% % Determine Decoders
% PHI = inv(Gamma)*Upsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #2: Using Element-Wise Operations - e.g., Matlab's "." syntax

% n = size(A,2); % number of neurons/decoders
% 
% %%% Determine Upsilon
% Upsilon = zeros(n,1);
% for ntor = 1:n
%     Upsilon(ntor) = sum(X.*A(:,ntor));
% end
% 
% %%% Determine Gamma
% Gamma = zeros(n,n);
% for itor = 1:n
%     for jtor = 1:n
%         Gamma(itor,jtor) = sum(A(:,itor).*A(:,jtor));
%     end
% end
% 
% % Determine Decoders
% PHI = inv(Gamma)*Upsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option #3: Using Vector Operations - e.g., the inner product

% Determine Upsilon
% - "similarity" between physical values and neural activation functions
Upsilon = A'*X;

% Determine Gamma
% - "similarity" among neural activation functions
Gamma = A'*A + 0.01*eye(size(A,2));


% Determine Decoders
PHI = Gamma\Upsilon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that there are at least two alternative methods for solving for PHI,
% and applying Eliasmith's equation ?.?. Some alternatives are:

% Determine Decoders - Alternative #1:
% PHI = (Gamma^-1)*Upsilon;

% Determine Decoders - Alternative #2:
% PHI = Gamma\Upsilon;

return
%eof