function [Cl, Circulacion, Circulacion_gradiente, Circulacion_estela] = CalculoNoEstacionario(Np, Ne, A, v, xpcontrol, xcgamma_estela, inc_t, xcgamma)
% CalculoNoEstacionario - Perform non-steady aerodynamic calculations
%
% Syntax: [Cl, Circulacion, Circulacion_gradiente, Circulacion_estela] =
%          CalculoNoEstacionario(Np, Ne, A, v, xpcontrol, xcgamma_estela, inc_t, xcgamma)
%
% Inputs:
%   Np - Number of profile divisions
%   Ne - Number of time steps
%   A - Coefficient matrix A
%   v - Input vector v
%   xpcontrol - Control point x-coordinate
%   xcgamma_estela - Estela point x-coordinate
%   inc_t - Time increment
%   xcgamma - x-coordinate vector
%
% Outputs:
%   Cl - Lift coefficient array
%   Circulacion - Circulation array
%   Circulacion_gradiente - Gradient of circulation array
%   Circulacion_estela - Circulation at estela point array
%
% Description: This function performs non-steady aerodynamic calculations
% using given inputs and coefficients. It calculates lift coefficients,
% circulation, gradient of circulation, and circulation at estela point.

% Initialize variables
RHS = ones(Np, 1);
A_estela = zeros(Np, Ne);
Circulacion_estela = zeros(1, Ne);
Circulacion = zeros(Np, Ne);
Circulacion_gradiente = zeros(Np, Ne);

% Loop over time steps
for i = 1:Ne
    for j = 1:Np
        A_estela(j, i) = -1 / (2 * pi) * (1 / (xpcontrol(j) - xcgamma_estela(i)));
    end
    if i == 1
        Matriz_A = inv(A - A_estela(:, 1) * transpose(RHS));
        Circulacion(:, i) = Matriz_A * v(:, i);
        Circulacion_estela(i) = -transpose(RHS) * Circulacion(:, i);
    else
        Aux = zeros(Np, i - 1);
        for k = 1:i - 1
            Aux(:, k) = (A_estela(:, 1) - A_estela(:, i - k + 1)) * Circulacion_estela(k);
        end
        Circulacion(:, i) = Matriz_A * (v(:, i) + sum(Aux, 2));
        Circulacion_estela(i) = -sum(Circulacion_estela) - transpose(RHS) * Circulacion(:, i);
    end
end

% Calculate circulation gradient
for i = 1:Np
    Circulacion_gradiente(i, :) = gradient(Circulacion(i, :)) / inc_t;
end

% Calculate lift coefficients
Cl = zeros(1, Ne);
Aux2 = zeros(Np, Ne);
for i = 1:Ne
    for j = 1:Np
        Aux2(j, i) = Circulacion_gradiente(j, i) * (1 - xcgamma(1, j));
    end
    Cl(1, i) = sum(Circulacion(:, i)) + sum(Aux2(:, i));
end

end
