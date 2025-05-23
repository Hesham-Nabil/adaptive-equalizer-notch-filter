function [g, y, E_cm] = cm_blind_equalizer(xn, M, mu, num_iter)
% CM_BLIND_EQUALIZER - Implements Constant Modulus Algorithm (CMA) Blind Equalization
% Inputs:
%   xn        - received signal (1D array)
%   M         - equalizer filter order
%   mu        - step size
%   num_iter  - number of iterations for adaptation
% Outputs:
%   g         - learned equalizer coefficients
%   y         - output of the equalizer over iterations
%   E_cm      - CM error per iteration

    g = ones(M+1, 1);            % Initial equalizer weights
    y = zeros(num_iter, 1);      % Equalizer output
    E_cm = zeros(num_iter, 1);   % CM error per iteration

    for n = M+1:num_iter
        x_vec = flipud(xn(n-M:n));     % Flipped input vector
        y(n) = g' * x_vec;             % Equalizer output
        E_cm(n) = (abs(y(n))^2 - 1)^2; % Constant modulus error
        g = g - mu * (abs(y(n))^2 - 1) * y(n) * x_vec; % Gradient descent
    end
end
