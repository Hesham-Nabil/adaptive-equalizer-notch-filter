function [a, y] = adaptive_filter(signal, mu, r)
% adaptive_filter - Adaptive Notch Filter using LMS
%
% Inputs:
%   signal : input signal vector
%   mu     : step size for LMS update
%   r      : pole radius (sharpness of the notch)
%
% Outputs:
%   a : vector of adaptive coefficient 'a[n]' over time
%   y : output filtered signal

% Initialization
N = length(signal);
a = zeros(N+1, 1);   % +1 to avoid overflow
y = zeros(N, 1);
e = zeros(N, 1);

% Adaptive filtering loop
for n = 1:N
    xn = signal(n);
    if n == 1
        xn_1 = 0; xn_2 = 0; yn_1 = 0; yn_2 = 0;
    elseif n == 2
        xn_1 = signal(n-1);
        xn_2 = 0;
        yn_1 = y(n-1);
        yn_2 = 0;
    else
        xn_1 = signal(n-1);
        xn_2 = signal(n-2);
        yn_1 = y(n-1);
        yn_2 = y(n-2);
    end

    % Error signal
    e(n) = xn + a(n) * xn_1 + xn_2;

    % Output signal
    y(n) = e(n) - r * a(n) * yn_1 - (r^2) * yn_2;

    % LMS update of 'a'
    a(n+1) = a(n) - mu * y(n) * xn_1;

    % Stability constraint
    if (a(n+1) > 2) || (a(n+1) < -2)
        a(n+1) = 0;
    end
end

% Remove last element to match signal length
a = a(1:N);

end
