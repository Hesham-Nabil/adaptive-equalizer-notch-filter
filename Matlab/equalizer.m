function [h, e, sn_delayed, y] = equalizer(sn, xn, hn, M, h_init)
% EQUALIZER Adaptive FIR equalizer using LMS algorithm
% Inputs:
%   sn: transmitted sequence (±1)
%   xn: received sequence (through channel + noise)
%   hn: channel impulse response
%   M: order of adaptive equalizer
%   h_init (optional): initial guess for the equalizer coefficients
%
% Outputs:
%   h: learned adaptive filter coefficients
%   e: error signal over time
%   sn_delayed: delayed version of input for comparison
%   y: equalizer output

    if nargin < 5
        h = ones(M+1,1);  % Default initialization
    else
        h = h_init;       % Custom initialization
    end

    u = 0.01;    % LMS step size
    N = length(sn);   % length of training sequence
    delays = round((length(hn)-1 + M)/2,0); % Approximate total system delay
    delays = 8; % override manually if you prefer (you did this in your code)
    
    hd = [zeros(1,delays), 1];
    sn_delayed = conv(hd,sn);
    
    e = zeros(length(sn),1);
    y = zeros(length(sn),1);
    
    % LMS training loop
    for n = 1:length(xn)
        xf  = flip(xn(1:n));
        % pad with zeros
        if length(xf) < M+1
            xf = [xf; zeros(M+1-length(xf),1)];
        end
        xf = xf(1:M+1);
        y(n) = sum(xf .* h);  % Equalizer output
        e(n) = sn_delayed(n) - y(n); % Error signal
        h = h + 2 * u * e(n) .* xf;  % Update the equalizer taps
    end
    
end
