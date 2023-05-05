% *Evaluating using Alpha Power and Theta Power*

addpath 'functions';

%% LOAD SIGNAL

data = csvread("my_matrix.csv");
Fs = 256;

vector_tp9 = data(1,:);
vector_af7 = data(2,:);
vector_af8 = data(3,:);
vector_tp10 = data(4,:);
vector_raux = data(5,:);

%% SIGNAL PRE-PROCESSING

vector_tp9 = vector_tp9 - vector_raux;
vector_af7 = vector_af7 - vector_raux;
vector_af8 = vector_af8 - vector_raux;
vector_tp10 = vector_tp10 - vector_raux;

low_cf = 1;
high_cf = 35;
filter_order = 4;
[b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
vector_tp9 = filtfilt(b,a,vector_tp9);
vector_af7 = filtfilt(b,a,vector_af7);
vector_af8 = filtfilt(b,a,vector_af8);
vector_tp10 = filtfilt(b,a,vector_tp10);

T = 1/Fs; % Sampling period
L = length(vector_tp9); % Length of signal
t = (0:L-1)*T; % Time vector
Y = fft(vector_tp9); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
window_size = 50;
P1 = movmean(P1, window_size);
f = Fs*(0:(L/2))/L; % Frequency vector
figure
plot(f,P1)
title('EEG FFT')

%% SIGNAL FILTERING

% Channel to be used
eeg_signal = [vector_tp9; vector_af8; vector_tp10];

% Define frequency bands of interest
alpha_freq = [8, 13];
theta_freq = [4, 8];
% Apply bandpass filter to signal
alpha_signal = eegfilt(eeg_signal, Fs, alpha_freq(1), alpha_freq(2));
theta_signal = eegfilt(eeg_signal, Fs, theta_freq(1), theta_freq(2));
%% ALPHA POWER SPECTRAL DENSITY (PSD)

window = 2*Fs;
alpha_power = zeros(size(eeg_signal,1), length(eeg_signal)-window+1);

for k = 1:size(alpha_signal,1)
    channel_signal = alpha_signal(k,:);
    for i = 1:length(alpha_power)
        singal_window = channel_signal(i:i+window-1);
    
        % Compute power spectral density
        [psd, freq] = pwelch(singal_window, [], [], [], Fs);
        
        % Integrate PSD over alpha frequency range to obtain alpha power
        alpha_freq_indices = freq >= alpha_freq(1) & freq <= alpha_freq(2);
        alpha_power(k, i) = trapz(freq(alpha_freq_indices), psd(alpha_freq_indices));
    end
end

figure
t = linspace(0, length(eeg_signal)/Fs, length(alpha_power(1,:)));
plot(t, alpha_power')
title('Alpha Power')

%% THETA POWER SPECTRAL DENSITY (PSD)

window = 2*Fs;
theta_power = zeros(size(eeg_signal,1), length(eeg_signal)-window+1);

for k = 1:size(theta_signal,1)
    channel_signal = theta_signal(k,:);
    for i = 1:length(alpha_power)
        singal_window = channel_signal(i:i+window-1);
    
        % Compute power spectral density
        [psd, freq] = pwelch(singal_window, [], [], [], Fs);
        
        % Integrate PSD over alpha frequency range to obtain alpha power
        theta_freq_indices = freq >= theta_freq(1) & freq <= theta_freq(2);
        theta_power(k, i) = trapz(freq(theta_freq_indices), psd(theta_freq_indices));
    end
end

figure
t = linspace(0, length(eeg_signal)/Fs,  length(theta_power(1,:)));
plot(t, theta_power')
title('Theta Power')