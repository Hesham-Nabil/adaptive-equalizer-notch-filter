clc; clear all; close all;

%% Part A1  :
%
% %Parameters
fs = 200;                  % Sampling frequency (Hz)
f_signal = 3;              % Desired signal frequency (Hz)
N = 500000;                % Number of samples
n = 0:N-1;                 % Time index
N_fft = 4096;              % FFT size

% Define experiments
experiments = [
    struct('mu',1e-6,'r',0.85,'f_noise',30),   % Experiment 1
    struct('mu',1e-5,'r',0.85,'f_noise',30),   % Experiment 2
    struct('mu',1e-6,'r',0.95,'f_noise',30),   % Experiment 3
    struct('mu',1e-6,'r',0.85,'f_noise',40)    % Experiment 4
    ];

legend_entries = {
    '\mu=1e-6, r=0.85, f=30Hz',
    '\mu=1e-5, r=0.85, f=30Hz',
    '\mu=1e-6, r=0.95, f=30Hz',
    '\mu=1e-6, r=0.85, f=40Hz'
    };

colors = {'b', 'r', 'g', 'm'};

% Initialize storage
X_f = cell(1,4);
Y_f = cell(1,4);
a_all = cell(1,4);
omega_est_all = cell(1,4);
omega_true_all = zeros(1,4);
y_all = cell(1,4);
x_all = cell(1,4);

% Loop through experiments
for idx = 1:length(experiments)
    mu = experiments(idx).mu;
    r = experiments(idx).r;
    f_noise = experiments(idx).f_noise;

    desired_signal = 1 * sin(2*pi*f_signal/fs * n);
    interference = 10 * sin(2*pi*f_noise/fs * n);
    x = desired_signal + interference;

    % Adaptive filtering
    a = zeros(N+1,1);
    y = zeros(N,1);
    e = zeros(N,1);

    for k = 1:N
        if k == 1
            x_n1 = 0; x_n2 = 0; y_n1 = 0; y_n2 = 0;
        elseif k == 2
            x_n1 = x(k-1); x_n2 = 0; y_n1 = y(k-1); y_n2 = 0;
        else
            x_n1 = x(k-1); x_n2 = x(k-2); y_n1 = y(k-1); y_n2 = y(k-2);
        end

        e(k) = x(k) + a(k) * x_n1 + x_n2;
        y(k) = e(k) - r * a(k) * y_n1 - r^2 * y_n2;
        a(k+1) = a(k) - mu * y(k) * x_n1;

        if (a(k+1) > 2) || (a(k+1) < -2)
            a(k+1) = 0;
        end
    end
    a = a(1:N);

    % Save results
    x_all{idx} = x;
    y_all{idx} = y;
    X_f{idx} = fftshift(fft(x(end-4096+1:end), N_fft));
    Y_f{idx} = fftshift(fft(y(end-4096+1:end), N_fft));
    a_all{idx} = a;
    omega_est_all{idx} = acos(-a/2);
    omega_true_all(idx) = 2*pi*f_noise/fs;
end

w = linspace(-pi, pi, N_fft);

% Plotting
for idx = 1:4
    figure;

    % Global title for the full figure
    sgtitle(['\mu=', num2str(experiments(idx).mu), ', r=', num2str(experiments(idx).r), ', f=', num2str(experiments(idx).f_noise), 'Hz']);

    % 1. Input vs Output Signal
    subplot(2,2,1); hold on;
    plot(n(end-50:end), x_all{idx}(end-50:end), '-', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'Input');
    plot(n(end-50:end), y_all{idx}(end-50:end), '--', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'Output');
    xlabel('n'); ylabel('Amplitude');
    title('Input and Output Signal');
    legend('show', 'FontSize', 8, 'Location', 'best');
    grid on;

    % 2. Magnitude Spectrum
    subplot(2,2,2); hold on;
    plot(w/pi, 20*log10(abs(X_f{idx})+1e-8), '-', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'Input');
    plot(w/pi, 20*log10(abs(Y_f{idx})+1e-8), '--', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'Output');
    xlabel('Normalized Frequency (\times \pi rad/sample)');
    ylabel('Magnitude (dB)');
    title('Magnitude Response of Input and Output');
    ylim([-100 100]);
    legend('show', 'FontSize', 8, 'Location', 'best');
    grid on;

    % 3. Evolution of a[n]
    subplot(2,2,3); hold on;
    plot(n, a_all{idx}, '-', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'a');
    plot(n, -2*cos(2*pi*experiments(idx).f_noise/fs)*ones(size(n)), '--', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'true a');
    xlabel('n'); ylabel('a[n]');
    title('Updating Parameter a');
    legend('show', 'FontSize', 8, 'Location', 'best');
    grid on;

    % 4. Tracking the Notch Frequency
    subplot(2,2,4); hold on;
    plot(n, omega_est_all{idx}, '-', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'w');
    plot(n, omega_true_all(idx)*ones(size(n)), '--', 'Color', colors{idx}, 'LineWidth', 1.2, 'DisplayName', 'true w');
    xlabel('n'); ylabel('\omega (rad/sample)');
    title('Updating the Notch Frequency');
    legend('show', 'FontSize', 8, 'Location', 'best');
    grid on;
end


%% Part A2  :
clc; clear all; close all;

% Parameters
fs = 200;                  % Sampling frequency (Hz)
f_signal = 3;              % Desired signal frequency (Hz)
N = 500000;                % Number of samples
n = 0:N-1;                 % Time index
N_fft = 4096;              % FFT size
mu = 1e-6;                 % Step size
r = 0.9;                   % Pole radius

%%Define slowly changing interference frequency
f_start = 30;              % Starting interference frequency (Hz)
f_end = 40;                % Ending interference frequency (Hz)
f_noise = linspace(f_start, f_end, N); % Linearly increasing frequency

% Desired signal
desired_signal = 1 * sin(2*pi*f_signal/fs * n);

% Interference signal with varying frequency
phi = cumsum(2*pi*f_noise/fs); % Instantaneous phase
interference = 10 * cos(phi);

% Combined input signal
x = desired_signal + interference;

% Adaptive filtering initialization
a = zeros(N+1,1);
y = zeros(N,1);
e = zeros(N,1);

for k = 1:N
    if k == 1
        x_n1 = 0; x_n2 = 0; y_n1 = 0; y_n2 = 0;
    elseif k == 2
        x_n1 = x(k-1); x_n2 = 0; y_n1 = y(k-1); y_n2 = 0;
    else
        x_n1 = x(k-1); x_n2 = x(k-2); y_n1 = y(k-1); y_n2 = y(k-2);
    end

    e(k) = x(k) + a(k)*x_n1 + x_n2;
    y(k) = e(k) - r*a(k)*y_n1 - r^2*y_n2;
    a(k+1) = a(k) - mu * y(k) * x_n1;

    if (a(k+1) > 2) || (a(k+1) < -2)
        a(k+1) = 0;
    end
end
a = a(1:N);

% Estimate omega
omega_est = acos(-a/2);
omega_true = 2*pi*f_noise/fs;

% plotting
figure;

% 1. Input and Output signals
subplot(2,2,1); hold on;
plot(n(end-50:end), x(end-50:end), '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Input');
plot(n(end-50:end), y(end-50:end), '--', 'Color', 'r', 'LineWidth', 1.8, 'DisplayName', 'Output');
xlabel('n'); ylabel('Amplitude');
title('Input and Output Signals (Last 50 samples)');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 2. Magnitude Spectrum
subplot(2,2,2); hold on;
plot(linspace(-fs/2, fs/2, N_fft), 20*log10(abs(fftshift(fft(x(end-4096+1:end),N_fft)))+1e-8), '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Input');
plot(linspace(-fs/2, fs/2, N_fft), 20*log10(abs(fftshift(fft(y(end-4096+1:end),N_fft)))+1e-8), '--', 'Color', 'r', 'LineWidth', 1.8, 'DisplayName', 'Output');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of Input and Output');
ylim([-100 100]);
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 3. Evolution of a[n]
subplot(2,2,3); hold on;
plot(n, a, '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'a');
plot(n, -2*cos(2*pi*f_noise/fs), '--', 'Color', 'r', 'LineWidth', 1.8, 'DisplayName', 'true a');
xlabel('n'); ylabel('a[n]');
title('Evolution of a[n]');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 4. Tracking Notch Frequency
subplot(2,2,4); hold on;
plot(n, omega_est, '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Estimated \\omega');
plot(n, omega_true, '--', 'Color', 'r', 'LineWidth', 1.8, 'DisplayName', 'True \\omega');
xlabel('n'); ylabel('\\omega (rad/sample)');
title('Tracking Notch Frequency');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% Optional: Add global title
sgtitle('Adaptive Notch Filter Tracking a Slowly Changing Frequency');




%% Part A3
clc; clear all; close all;

% Parameters
fs = 200;                  % Sampling frequency (Hz)
N = 50000;                 % Number of samples
n = 0:N;                   % Time vector
u = 1e-5;                  % Step size
r = 0.95;                  % Pole radius

% Signal definition
fsig = 5;                  % Desired signal frequency (Hz)
fn1 = 45;                  % Interference frequency 1 (Hz)
fn2 = 60;                  % Interference frequency 2 (Hz)

desired = sin(2*pi*fsig/fs * n);         % Desired signal
noise1 = 8*sin(2*pi*fn1/fs * n);          % Interference 1
noise2 = 4*sin(2*pi*fn2/fs * n);          % Interference 2
x = desired + noise1 + noise2;            % Observed input signal

% True a values
wnoise1 = 2*pi*fn1/fs;
wnoise2 = 2*pi*fn2/fs;
atrue1 = -2*cos(wnoise1);
atrue2 = -2*cos(wnoise2);

% Initialization
a1 = zeros(length(x),1);
a2 = zeros(length(x),1);
yn1 = zeros(length(x),1);
yn2 = zeros(length(x),1);

% Adaptive Filtering
for k = 1:length(x)
    [yn1, xn1_1] = adaptive_filter_OH(x, a1, yn1, u, r, k);
    [yn2, xn2_1] = adaptive_filter_OH(yn1, a2, yn2, u, r, k);

    % Update a1
    a1(k+1) = a1(k) - u * yn1(k) * xn1_1;
    if (a1(k+1) > 2) || (a1(k+1) < -2)
        a1(k+1) = 0;
    end

    % Update a2
    a2(k+1) = a2(k) - u * yn2(k) * xn2_1;
    if (a2(k+1) > 2) || (a2(k+1) < -2)
        a2(k+1) = 0;
    end
end

a1 = a1(1:end-1);
a2 = a2(1:end-1);

% Frequency Analysis
N_fft = 512;
w = linspace(-pi, pi, N_fft);
X_f = fftshift(fft(x, N_fft));
Y1_f = fftshift(fft(yn1(end-511:end), N_fft));
Y2_f = fftshift(fft(yn2(end-511:end), N_fft));

% Plotting
figure;

% 1. Magnitude spectrum of x[n]
subplot(3,2,1);
plot(w, abs(X_f), 'b');
title('Magnitude Response of x[n]');
xlabel('\omega (rad/sample)');
ylabel('Magnitude');
grid on;

% 2. Magnitude spectrum of y1[n]
subplot(3,2,3);
plot(w, abs(Y1_f), 'b');
title('Magnitude Response of y1[n]');
xlabel('\omega (rad/sample)');
ylabel('Magnitude');
grid on;

% 3. Magnitude spectrum of y2[n]
subplot(3,2,5);
plot(w, abs(Y2_f), 'b');
title('Magnitude Response of y2[n]');
xlabel('\omega (rad/sample)');
ylabel('Magnitude');
grid on;

% 4. x[n] vs y1[n] vs y2[n]
subplot(3,2,6); hold on;
plot(n(end-30:end), x(end-30:end), 'b-', 'LineWidth', 1);
plot(n(end-30:end), yn1(end-30:end), 'r--', 'LineWidth', 1);
plot(n(end-30:end), yn2(end-30:end), 'm-', 'LineWidth', 1);
legend('x[n]', 'y1[n]', 'y2[n]');
title('Comparison of x[n], y1[n], and y2[n]');
xlabel('n');
ylabel('Amplitude');
grid on;
ylim([-20 20]);

% 5. Evolution of a1[n] and true a1
subplot(3,2,2); hold on;
plot(a1, 'b-', 'LineWidth', 1);
plot(atrue1*ones(length(a1),1), 'r--', 'LineWidth', 1);
legend('a1', 'true a1');
title('Updating Parameter a1');
xlabel('n');
ylabel('a_1[n]');
grid on;

% 6. Evolution of a2[n] and true a2
subplot(3,2,4); hold on;
plot(a2, 'b-', 'LineWidth', 1);
plot(atrue2*ones(length(a2),1), 'r--', 'LineWidth', 1);
legend('a2', 'true a2');
title('Updating Parameter a2');
xlabel('n');
ylabel('a_2[n]');
grid on;

sgtitle('Adaptive Notch Filter for Two Sinusoidal Interferences');

% Colors
input_color = 'b'; % blue
output_color = 'r'; % red

figure;

% 1. Input and Output Signals (last 50 samples)
subplot(2,2,1); hold on;
plot(n(end-50:end), x(end-50:end), input_color, 'LineWidth', 1.5, 'DisplayName', 'Input');
plot(n(end-50:end), yn2(end-50:end), output_color, 'LineWidth', 1.5, 'DisplayName', 'Output');
xlabel('n'); ylabel('Amplitude');
title('Input and Output Signals (Last 50 samples)');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 2. Magnitude Spectrum
subplot(2,2,2); hold on;
N_fft = 4096;
plot(linspace(-fs/2, fs/2, N_fft), 20*log10(abs(fftshift(fft(x(end-4096+1:end),N_fft)))+1e-8), input_color, 'LineWidth', 1.5, 'DisplayName', 'Input');
plot(linspace(-fs/2, fs/2, N_fft), 20*log10(abs(fftshift(fft(yn2(end-4096+1:end),N_fft)))+1e-8), output_color, 'LineWidth', 1.5, 'DisplayName', 'Output');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Spectrum of Input and Output');
ylim([-100 100]);
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 3. Evolution of a1[n] and a2[n]
subplot(2,2,3); hold on;
plot(n, a1, input_color, 'LineWidth', 1.5, 'DisplayName', 'a_1');
plot(n, a2, output_color, 'LineWidth', 1.5, 'DisplayName', 'a_2');
plot(n, ones(size(n))*atrue1, input_color, 'LineWidth', 1, 'DisplayName', 'true a_1');
plot(n, ones(size(n))*atrue2, output_color, 'LineWidth', 1, 'DisplayName', 'true a_2');
xlabel('n');
ylabel('a[n]');
title('Evolution of a_1[n] and a_2[n]');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% 4. Tracking Notch Frequencies
subplot(2,2,4); hold on;
omega1_est = acos(-a1/2);
omega2_est = acos(-a2/2);
omega1_true = 2*pi*fn1/fs;
omega2_true = 2*pi*fn2/fs;


plot(n, omega1_est, '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Estimated \omega_1');
plot(n, omega2_est, '--', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Estimated \omega_2');
plot(n, omega1_true*ones(size(n)), input_color, 'LineWidth', 1, 'DisplayName', 'True \omega_1');
plot(n, omega2_true*ones(size(n)), output_color, 'LineWidth', 1, 'DisplayName', 'True \omega_2');
xlabel('n');
ylabel('\omega (rad/sample)');
title('Tracking Notch Frequencies');
legend('show', 'FontSize',8, 'Location', 'best');
grid on;

% Global title
sgtitle('Adaptive Notch Filter for Two Interferences');



%% Part B: Adaptive Equalization 
    clc; clear all; close all;
    
    % Generate Training Data
    rng(1, 'twister');
    sn = randi([0, 1], 1, 1000);
    sn(sn == 0) = -1;
    
    hn = [0.3, 1, 0.7, 0.3, 0.2]; % Channel impulse response
    noise = wgn(1, length(sn) + length(hn) - 1, -25); % Additive Gaussian noise
    xn = (conv(sn, hn) + noise)';
    
    % Parameters
    M = 11; % Equalizer order
    delays = 8; % Fixed delay
    h = ones(M+1, 1); % Initial filter coefficients
    u = 0.01; % Step size
    
    hd = [zeros(1, delays), 1]; % Delay line
    sn_delayed = conv(hd, sn)';
    
    e = zeros(length(sn), 1);
    y = zeros(length(sn), 1);
    h_track = zeros(length(sn), M+1);
    overall_lti_track = zeros(length(sn), M+5);
    
    % Adaptive LMS Training
    for n = 1:length(xn)
        xf = flip(xn(max(1, n - M):n));
        xf = [xf; zeros(M+1-length(xf), 1)];
        y(n) = sum(xf .* h);
        e(n) = sn_delayed(n) - y(n);
        h = h + 2 * u * e(n) .* xf;
        h_track(n, :) = h';
        overall_lti_track(n, :) = conv(h, hn);
    end
    
    overall_lti = conv(h, hn);
    
    % Training Error and Output Comparison
    figure;
    subplot(2,1,1);
    plot(abs(e), 'Color', [0.1, 0.6, 0.8]);
    xlabel('n'); ylabel('|Error|');
    title('Absolute Value of Error during Training');
    
    subplot(2,1,2);
    stem(900:1004, sn_delayed(900:1004), 'filled', 'MarkerFaceColor', 'b');
    hold on;
    stem(900:1004, y(900:1004), 'filled', 'MarkerFaceColor', 'm');
    legend('Delayed Input', 'Equalizer Output');
    xlabel('n'); ylabel('Signal');
    title('Comparison: Input vs Output after Equalization');
    saveas(gcf, 'figures/training_error_comparison.png');
    
    % Accuracy Calculation
    in = sn_delayed(900:1004);
    out = y(900:1004);
    out(out > 0) = 1;
    out(out < 0) = -1;
    out(out == 0) = -1;
    acc = (1 - abs(sum(in' - out)) / length(in)) * 100;
    sprintf('Accuracy = %.2f%%', acc)
    
    % Filter Convergence Plots
    figure;
    for i = 1:M+1
        subplot(3,4,i);
        plot(h_track(:,i), 'Color', [rand, rand, rand]);
        xlabel('n');
        title(sprintf('Tap %d', i-1));
    end
    sgtitle('Convergence of Equalizer Coefficients');
    saveas(gcf, 'figures/equalizer_convergence.png');
    
    figure;
    for i = 1:M+5
        subplot(4,4,i);
        plot(overall_lti_track(:,i), 'Color', [rand, rand, rand]);
        hold on;
        if i ~= 9
            plot(zeros(1,length(sn)),'--r');
        else
            plot(ones(1,length(sn)),'--g');
        end
        xlabel('n');
        title(sprintf('Tap %d', i-1));
    end
    sgtitle('Convergence of Overall LTI System');
    saveas(gcf, 'figures/overall_lti_convergence.png');
    
    % Frequency Response and Impulse Response
    figure;
    subplot(3,2,1);
    [channel_f, w] = DTFT(hn);
    plot(w, abs(channel_f), 'LineWidth', 1.5);
    xlabel('\omega (rad/s)'); ylabel('Magnitude');
    title('Channel Frequency Response');
    
    subplot(3,2,2);
    plot(w, angle(channel_f), 'LineWidth', 1.5);
    xlabel('\omega (rad/s)'); ylabel('Phase');
    title('Channel Phase Response');
    
    subplot(3,2,3);
    [h_f, ~] = DTFT(h);
    plot(w, abs(h_f), 'LineWidth', 1.5);
    xlabel('\omega (rad/s)'); ylabel('Magnitude');
    title('Equalizer Frequency Response');
    
    subplot(3,2,4);
    plot(w, angle(h_f), 'LineWidth', 1.5);
    xlabel('\omega (rad/s)'); ylabel('Phase');
    title('Equalizer Phase Response');
    
    subplot(3,2,5);
    plot(w, abs(h_f .* channel_f'), 'LineWidth', 1.5);
    hold on;
    plot(w, ones(1,length(w)), '--r');
    xlabel('\omega (rad/s)'); ylabel('Magnitude');
    title('Overall System Frequency Response');
    
    subplot(3,2,6);
    plot(w, angle(h_f .* channel_f'), 'LineWidth', 1.5);
    xlabel('\omega (rad/s)'); ylabel('Phase');
    title('Overall System Phase Response');
    saveas(gcf, 'figures/frequency_response.png');
    
    figure;
    subplot(3,1,1);
    stem(0:length(hn)-1, hn, 'filled');
    xlabel('n'); ylabel('Amplitude');
    title('Channel Impulse Response');
    
    subplot(3,1,2);
    stem(0:length(h)-1, h, 'filled');
    xlabel('n'); ylabel('Amplitude');
    title('Equalizer Impulse Response');
    
    subplot(3,1,3);
    stem(0:length(overall_lti)-1, overall_lti, 'filled');
    xlabel('n'); ylabel('Amplitude');
    title('Overall LTI System Impulse Response');
    saveas(gcf, 'figures/impulse_responses.png');
    
    figure;
    subplot(3,1,1);
    zplane(hn, [1]);
    title('Pole-Zero plot: Channel');
    
    subplot(3,1,2);
    zplane(h, [1]);
    title('Pole-Zero plot: Equalizer');
    
    subplot(3,1,3);
    zplane(conv(h, hn), [1]);
    title('Pole-Zero plot: Overall System');
    saveas(gcf, 'figures/pole_zero_plots.png');
    
    % Testing the Equalizer (new input)
    rng(5, 'twister');
    sn = randi([0,1],1,1000);
    sn(sn==0) = -1;
    noise = wgn(1,length(sn)+length(hn)-1,-25);
    xn = (conv(sn,hn)+noise)';
    yn = conv(h,xn);
    
    % Hard Decision
    yn(yn>0) = 1; yn(yn<0) = -1; yn(yn==0) = -1;
    
    % For comparison
    hd = [zeros(1,8),1];
    sn_delayed = conv(hd, sn)';
    
    % Plot Results
    figure;
    subplot(2,2,1);
    stem(908:1008, sn_delayed(908:1008), 'filled', 'MarkerFaceColor', 'k');
    hold on;
    stem(908:1008, yn(908:1008), 'filled', 'MarkerFaceColor', 'c');
    legend('Input','Equalized Output');
    xlabel('n'); ylabel('Signal'); title('Equalized Input vs Output');
    
    subplot(2,2,2);
    err = sn_delayed(908:1008) - yn(908:1008);
    stem(908:1008, err, 'MarkerFaceColor', 'r');
    c = sum(err == 0) / length(err) * 100;
    title(sprintf('Error (Accuracy = %.2f%%)', c));
    xlabel('n'); ylabel('Error');
    
    subplot(2,2,3);
    sn_delayed2 = conv([0 1], sn)';
    xn(xn>0) = 1; xn(xn<0) = -1; xn(xn==0) = -1;
    stem(901:1001, sn_delayed2(901:1001), 'filled', 'MarkerFaceColor', 'g');
    hold on;
    stem(901:1001, xn(901:1001), 'filled', 'MarkerFaceColor', 'm');
    legend('Input','Un-equalized Output');
    xlabel('n'); ylabel('Signal'); title('Un-equalized Input vs Output');
    
    subplot(2,2,4);
    err = sn_delayed2(901:1001) - xn(901:1001);
    stem(901:1001, err, 'MarkerFaceColor', 'b');
    c = sum(err == 0) / length(err) * 100;
    title(sprintf('Error (Accuracy = %.2f%%)', c));
    xlabel('n'); ylabel('Error');
    saveas(gcf, 'figures/test_equalizer_performance.png');
    
    
    
    %Performance Analysis - SNR, Filter Order, Initialization
    
    clc; clear all; close all;
    
    % Generate Training Data
    rng(1,'twister');
    sn = randi([0,1], 1,1000);
    sn(sn==0) = -1;
    hn = [0.3, 1, 0.7, 0.3, 0.2];
    
    noise = wgn(1,length(sn)+length(hn)-1,-25);
    xn = (conv(sn,hn) + noise)';
    [h,e,sn_delayed,out] = equalizer(sn,xn,hn,11);
    
    % Test Setup
    rng(5,'twister');
    sn = randi([0,1], 1,1000);
    sn(sn==0) = -1;
    delays = 8;
    hd = [zeros(1,delays), 1];
    sn_delayed = conv(hd,sn);
    
    % Performance vs SNR
    err1 = zeros(61,1);
    acc = zeros(61,1);
    for i = 0:60
        noise = wgn(1,length(sn)+length(hn)-1,-i);
        xn = (conv(sn,hn) + noise)';
        out = conv(h,xn);
        outd = sign(out);
        outd(outd==0) = -1;
        err1(i+1) = sum(abs(sn_delayed(8:900)' - out(8:900)));
        acc(i+1) = (1 - abs(sum(sn_delayed(8:900)' - out(8:900))) / length(8:900)) * 100;
    end
    
    figure;
    subplot(2,1,1);
    plot(0:60, err1, 'Color', [0.1 0.6 0.8], 'LineWidth', 1.5);
    hold on;
    plot(0:60, ones(1,61)*167.94, '--r', 'LineWidth', 1.5);
    xlabel('SNR (dB)'); ylabel('Sum of |Error|');
    title('Equalization Error vs SNR');
    legend('Error', 'Error (No noise)');
    grid on;
    
    subplot(2,1,2);
    plot(0:60, acc, 'Color', [0.8 0.3 0.1], 'LineWidth', 1.5);
    hold on;
    plot(0:60, ones(1,61)*99.84, '--r', 'LineWidth', 1.5);
    xlabel('SNR (dB)'); ylabel('Accuracy (%)');
    title('Equalization Accuracy vs SNR');
    legend('Accuracy', 'Accuracy (No noise)');
    grid on;
    
    % Performance vs Filter Order
    clc; clear err1; clear acc;
    noise = wgn(1,length(sn)+length(hn)-1,-25);
    xn = (conv(sn,hn) + noise)';
    
    err1 = zeros(21,1);
    acc = zeros(21,1);
    for i = 0:20
        [h,e,sn_delayed,out] = equalizer(sn,xn,hn,i+3);
        out = conv(h,xn);
        outd = sign(out);
        outd(outd==0) = -1;
        err1(i+1) = sum(abs(sn_delayed(200:900)' - out(200:900)));
        acc(i+1) = (1 - abs(sum(sn_delayed(200:900)' - out(200:900))) / length(8:900)) * 100;
    end
    
    figure;
    subplot(2,1,1);
    plot(3:23, err1, 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    xlabel('Filter Order M'); ylabel('Sum of |Error|');
    title('Equalization Error vs Filter Order');
    grid on;
    
    subplot(2,1,2);
    plot(3:23, acc, 'Color', [0.6 0.2 0.8], 'LineWidth', 1.5);
    xlabel('Filter Order M'); ylabel('Accuracy (%)');
    title('Equalization Accuracy vs Filter Order');
    grid on;
    
    % Performance vs Initialization
    clc; clear err1; clear acc;
    noise = wgn(1,length(sn)+length(hn)-1,-25);
    xn = (conv(sn,hn) + noise)';
    
    c = 1;
    err1 = zeros(4001,1);
    acc = zeros(4001,1);
    for i = -2000:2000
        [h,e,sn_delayed,out] = equalizer(sn,xn,hn,11,ones(12,1)*i);
        out = conv(h,xn);
        outd = sign(out);
        outd(outd==0) = -1;
        err1(c) = sum(abs(sn_delayed(200:900)' - out(200:900)));
        acc(c) = (1 - abs(sum(sn_delayed(200:900)' - out(200:900))) / length(8:900)) * 100;
        c = c + 1;
    end
    
    figure;
    subplot(2,1,1);
    plot(-2000:2000, err1, 'Color', [0.4 0.4 0.9], 'LineWidth', 1.5);
    xlabel('Initial Tap Value'); ylabel('Sum of |Error|');
    title('Error vs Initial Equalizer Tap Value');
    grid on;
    
    subplot(2,1,2);
    plot(-2000:2000, acc, 'Color', [0.9 0.4 0.4], 'LineWidth', 1.5);
    xlabel('Initial Tap Value'); ylabel('Accuracy (%)');
    title('Accuracy vs Initial Equalizer Tap Value');
    grid on;

%% Part C
clc; clear all; close all;

% Input Signal
rng(1,'twister');
sn = randi([0,1], 1,100000);
sn(sn==0) = -1;

% Channel and Noise
h_channel = [0.3, 1, 0.7, 0.3, 0.2];
noise = wgn(1,length(sn)+length(h_channel)-1,-25);
xn = conv(sn,h_channel)';

% Initialization
M = 11;
h = ones(M+1,1);
u = 0.0002;
e = zeros(length(sn),1);
y = zeros(length(sn),1);
h_track = zeros(length(sn),M+1);
overall_lti_track = zeros(length(sn),M+5);

% CMA Equalizer Update
for n = 1:length(xn)
    xf = flip(xn(1:n));
    if length(xf) < M+1
        xf = [xf; zeros(M+1-length(xf),1)];
    end
    xf = xf(1:M+1);
    y(n) = sum(xf .* h);
    e(n) = y(n)^2 - 1;
    h = h - 4 * u * e(n) * y(n) .* xf;
    h_track(n,:) = h';
    overall_lti_track(n,:) = conv(h,h_channel);
end

overall_lti = conv(h,h_channel);

% Accuracy Evaluation
mstart = length(sn)-1000;
mend = length(sn)-900;
delays = M-4;
hd = [zeros(1,delays), 1];
sn_delayed = conv(hd,sn);
in = sn_delayed(mstart:mend);
out = sign(y(mstart:mend))';
a = in + out;
acc = sum(a==0) / length(a) * 100;
sprintf('Accuracy = %d%%', round(acc, 0))

% Error Plot
figure;
plot(abs(e));
xlabel('n'); ylabel('|Error|');
title('Absolute Value of the Error');
saveas(gcf, 'figures2/error_abs.png');

% Input vs Output Comparison
figure;
subplot(2,1,1);
stem(0:mend-mstart, sn_delayed(mstart:mend)); hold on;
stem(0:mend-mstart, y(mstart:mend));
legend('Input','Output');
title('Input (s[n]) vs Output (y[n])'); xlabel('n'); ylabel('Signal');

subplot(2,1,2);
stem(0:mend-mstart, sn_delayed(mstart:mend)); hold on;
stem(0:mend-mstart, sign(y(mstart:mend)));
legend('Input','Output');
title(sprintf('Input vs Output after Decision Block (%d%% accuracy)', round(acc, 0)));
xlabel('n'); ylabel('Signal');
saveas(gcf, 'figures2/input_vs_output.png');

% Output Amplitude
figure;
plot(abs(y)); hold on;
plot(ones(1,length(y)), '--r','LineWidth',1.5);
title('Output Amplitude'); xlabel('n');
saveas(gcf, 'figures2/output_amplitude.png');

% Equalization vs No Equalization
figure;
subplot(2,2,1);
delays = M-4;
hd = [zeros(1,delays), 1];
sn_delayed = -1 * conv(hd,sn);
stem(0:mend-mstart, sn_delayed(mstart:mend)); hold on;
stem(0:mend-mstart, sign(y(mstart:mend)));
legend('Input', 'Output (inverted)');
title('Input vs Output after Equalization'); xlabel('n'); ylabel('Signal');

subplot(2,2,2);
a = sn_delayed(mstart:mend) - sign(y(mstart:mend))';
stem(0:mend-mstart, a);
c = sum(a==0) / length(a) * 100;
title(sprintf('Error (accuracy = %d%%) with Equalization', round(c, 0)));
xlabel('n'); ylabel('Error');

subplot(2,2,3);
delays = 1;
hd = [zeros(1,delays), 1];
sn_delayed = conv(hd,sn);
stem(0:mend-mstart, sn_delayed(mstart:mend)); hold on;
stem(0:mend-mstart, sign(xn(mstart:mend)));
legend('Input','Output');
title('Input vs Output (Un-equalized)'); xlabel('n'); ylabel('Signal');

subplot(2,2,4);
stem(0:mend-mstart, sn_delayed(mstart:mend) - sign(xn(mstart:mend)));
c = sum(sn_delayed(mstart:mend) == sign(xn(mstart:mend))) / length(a) * 100;
title(sprintf('Error (accuracy = %d%%) Un-equalized', round(c, 0)));
xlabel('n'); ylabel('Error');
saveas(gcf, 'figures2/equalization_comparison.png');

% Convergence of Equalizer
figure;
for i = 1:M+1
    subplot(3,4,i);
    plot(h_track(:,i));
    xlabel('n'); title(sprintf('Tap %d', i-1));
end
sgtitle('Equalizer Convergence Over Time');
saveas(gcf, 'figures2/equalizer_convergence.png');

% Convergence of Overall System
figure;
for i = 1:M+5
    subplot(4,4,i);
    plot(overall_lti_track(:,i)); hold on;
    if i ~= 8
        plot(zeros(1,length(sn)), '--r');
    else
        plot(-1*ones(1,length(sn)), '--r');
    end
    xlabel('n'); title(sprintf('Tap %d', i-1));
end
sgtitle('Overall LTI System Convergence');
saveas(gcf, 'figures2/lti_convergence.png');

% Frequency and Impulse Response
figure;
subplot(3,2,1);
[channel_f, w] = DTFT(h_channel);
plot(w, abs(channel_f));
xlabel('\omega (rad/s)'); ylabel('Gain'); xlim([0, pi]);
title('Channel Frequency Magnitude');

subplot(3,2,2);
plot(w, angle(channel_f));
xlabel('\omega (rad/s)'); ylabel('Phase'); xlim([0, pi]);
title('Channel Frequency Phase');

subplot(3,2,3);
[h_f, ~] = DTFT(h);
plot(w, abs(h_f));
xlabel('\omega (rad/s)'); ylabel('Gain'); xlim([0, pi]);
title('Equalizer Frequency Magnitude');

subplot(3,2,4);
plot(w, angle(h_f));
xlabel('\omega (rad/s)'); ylabel('Phase'); xlim([0, pi]);
title('Equalizer Frequency Phase');

subplot(3,2,5);
plot(w, abs(h_f .* channel_f'));
hold on;
plot(w, ones(1,length(w)), '--r');
xlabel('\omega (rad/s)'); ylabel('Gain'); xlim([0, pi]);
title('Overall Frequency Magnitude');

subplot(3,2,6);
plot(w, angle(h_f .* channel_f'));
xlabel('\omega (rad/s)'); ylabel('Phase'); xlim([0, pi]);
title('Overall Frequency Phase');
saveas(gcf, 'figures2/frequency_response.png');

% Impulse Responses
figure;
subplot(3,1,1);
stem(0:length(h_channel)-1, h_channel);
xlabel('n'); ylabel('h[n]');
title('Channel Impulse Response');

subplot(3,1,2);
stem(0:length(h)-1, h);
xlabel('n'); ylabel('h[n]');
title('Equalizer Impulse Response');

subplot(3,1,3);
a_lti = conv(h,h_channel);
stem(0:length(a_lti)-1, a_lti);
xlabel('n'); ylabel('h[n]');
title('Overall LTI System Impulse Response');
saveas(gcf, 'figures2/impulse_response.png');
