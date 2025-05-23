%% Part B - Adaptive Equalization (Student Style, Corrected)
clc; clear all; close all;

% Parameters
N = 1000;                   % Number of symbols
M = 11;                     % Equalizer filter order
mu = 0.01;                  % LMS step size
h_channel = [0.3, 1, 0.7, 0.3, 0.2]; % Channel impulse response
delay = 8;                  % Fixed delay (student uses 8)

%% --- Training Phase ---

% Generate random input
rng(1, 'twister');
s = randi([0,1], 1, N);
s(s==0) = -1;

% Pass through channel + add noise
xn_clean = conv(s, h_channel);
noise = wgn(1,length(xn_clean),-25);  % -25 dB noise (very low)
xn = xn_clean + noise;                % Received noisy signal

% Generate delayed desired signal
hd = [zeros(1,delay), 1];
s_delayed = conv(hd, s);
s_delayed = s_delayed(1:length(xn));  % Match size

% Initialize adaptive filter
h_eq = ones(M+1,1);    % Initial equalizer taps
y = zeros(length(xn),1);
e = zeros(length(xn),1);
h_track = zeros(length(xn), M+1);
overall_lti_track = zeros(length(xn), M + length(h_channel));

% LMS Adaptation
for n = 1:length(xn)
    if n <= M
        xf = flip([xn(1:n)'; zeros(M+1-n,1)]);  % Pad then flip
    else
        xf = flip(xn(n-M:n)');
    end

    y(n) = h_eq' * xf;
    e(n) = s_delayed(n) - y(n);

    h_eq = h_eq + 2 * mu * e(n) * xf; % LMS Update
    h_track(n,:) = h_eq;
    overall_lti_track(n,:) = conv(h_eq', h_channel);
end

overall_lti = conv(h_eq', h_channel);

%% --- Plots for Training Phase ---

% 1. Absolute Error
figure;
subplot(2,1,1);
plot(abs(e));
xlabel('n'); ylabel('abs(error)');
title('Absolute value of error');

subplot(2,1,2);
stem(900:1004, s_delayed(900:1004), 'filled'); hold on;
stem(900:1004, y(900:1004), 'r');
legend('Delayed input','Equalizer output');
title('Input vs Output after training');
xlabel('n'); ylabel('Amplitude');

% 2. Convergence of Equalizer Coefficients
figure;
for k = 1:M+1
    subplot(3,4,k);
    plot(h_track(:,k));
    xlabel('n');
    title(sprintf('Tap %d',k-1));
end
sgtitle('Equalizer Coefficient Convergence');

% 3. Convergence of Overall System Coefficients
figure;
for k = 1:length(overall_lti)
    subplot(4,4,k);
    plot(overall_lti_track(:,k));
    hold on;
    if k == delay+1
        plot(ones(1,length(xn)),'r--');
    else
        plot(zeros(1,length(xn)),'r--');
    end
    xlabel('n');
    title(sprintf('Overall Tap %d',k-1));
end
sgtitle('Overall System Convergence');

% 4. Frequency Response (optional)
figure;
freqz(h_eq,1,512);
title('Frequency Response of Equalizer');

%% --- Testing Phase ---

% New random input
rng(5,'twister');
s_new = randi([0,1],1,N);
s_new(s_new==0) = -1;

% Pass through channel + noise
xn_new_clean = conv(s_new, h_channel);
noise_new = wgn(1,length(xn_new_clean),-25);
xn_new = xn_new_clean + noise_new;

% Desired delayed signal
s_new_delayed = conv(hd, s_new);
s_new_delayed = s_new_delayed(1:length(xn_new));

% Apply equalizer
yn = conv(h_eq', xn_new);
yn = yn(1:length(xn_new));

% Decision block
yn_decision = sign(yn);
xn_new_decision = sign(xn_new);

yn_decision(yn_decision==0) = -1;
xn_new_decision(xn_new_decision==0) = -1;

% Accuracy calculation
accuracy_eq = mean(yn_decision(900:1000)' == s_new_delayed(900:1000)') * 100;
accuracy_uneq = mean(xn_new_decision(900:1000)' == s_new_delayed(900:1000)') * 100;

fprintf('Equalized Accuracy: %.2f%%\n', accuracy_eq);
fprintf('Un-equalized Accuracy: %.2f%%\n', accuracy_uneq);

%% --- Testing Plots

figure;
subplot(2,2,1);
stem(900:1000, s_new_delayed(900:1000),'filled'); hold on;
stem(900:1000, yn_decision(900:1000),'r--');
legend('Input', 'Equalized Output');
title('Input vs Equalized Output');

subplot(2,2,2);
stem(900:1000, s_new_delayed(900:1000) - yn_decision(900:1000));
title(sprintf('Equalized Error (Accuracy=%.1f%%)', accuracy_eq));

subplot(2,2,3);
stem(900:1000, s_new_delayed(900:1000),'filled'); hold on;
stem(900:1000, xn_new_decision(900:1000),'r--');
legend('Input', 'Un-equalized Output');
title('Input vs Un-equalized Output');

subplot(2,2,4);
stem(900:1000, s_new_delayed(900:1000) - xn_new_decision(900:1000));
title(sprintf('Un-equalized Error (Accuracy=%.1f%%)', accuracy_uneq));

sgtitle('Equalizer Testing on New Input Sequence');

