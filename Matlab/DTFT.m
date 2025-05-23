function [X, w] = DTFT(x)
% DTFT - Compute the Discrete Time Fourier Transform of a signal
%
% Usage:
% [X, w] = DTFT(x)
%
% Inputs:
%   x : input discrete-time signal
%
% Outputs:
%   X : DTFT of the signal
%   w : corresponding frequency vector (radians/sample)

N = length(x);
N_fft = 4096; % you can adjust to higher for better resolution if needed
X = fft(x, N_fft);
X = fftshift(X); % shift zero frequency to center
w = linspace(-pi, pi, N_fft);
end
