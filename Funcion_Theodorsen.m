function [C] = Funcion_Theodorsen(k)
% Funcion_Theodorsen - Calculate Theodorsen's function for a given parameter k
%
% Syntax: [C] = Funcion_Theodorsen(k)
%
% Inputs:
%   k - Parameter value
%
% Outputs:
%   C - Theodorsen's function value
%
% Description: This function calculates Theodorsen's function using Bessel
%              functions and returns its complex value for a given parameter k.

% Calculate Bessel functions
H12 = besselh(1, 2, k);
H02 = besselh(0, 2, k);

% Calculate Theodorsen's function
C = complex(H12 / (H12 + (1i) * H02));

end
