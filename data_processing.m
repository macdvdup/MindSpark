clear all; close all;
addpath 'functions'

data = csvread("my_matrix.csv");
Fs = 256; % Sampling frequency
data = data(:,:);

vector_tp9 = data(1,:);
vector_af7 = data(2,:);
vector_af8 = data(3,:);
vector_tp10 = data(4,:);
vector_raux = data(5,:);

vector_tp9 = vector_tp9 - vector_raux;
vector_af7 = vector_af7 - vector_raux;
vector_af8 = vector_af8 - vector_raux;
vector_tp10 = vector_tp10 - vector_raux;

low_cf = 1;
high_cf = 35;
filter_order = 6;
[b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
vector_tp9 = filtfilt(b,a,vector_tp9);

plot(vector_tp9);

T = 1/Fs; % Sampling period
L = length(vector_tp9); % Length of signal
t = (0:L-1)*T; % Time vector
Y = fft(vector_tp9); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L; % Frequency vector
figure
plot(f,P1) 

% Load EEG data
eeg_data = vector_tp9;

% Define frequency bands of interest
alpha_freq = [8, 13];
theta_freq = [4, 8];
beta_freq = [13, 30];
gamma_freq = [30, 100];

% Apply bandpass filter to signal
alpha_signal = eegfilt(eeg_data, Fs, alpha_freq(1), alpha_freq(2));
theta_signal = eegfilt(eeg_data, Fs, theta_freq(1), theta_freq(2));
beta_signal = eegfilt(eeg_data, Fs, beta_freq(1), beta_freq(2));
gamma_signal = eegfilt(eeg_data, Fs, gamma_freq(1), gamma_freq(2));

figure;
subplot(4,1,1);
plot(alpha_signal);
title('Alpha band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,2);
plot(theta_signal);
title('Theta band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,3);
plot(beta_signal);
title('Beta band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,4);
plot(gamma_signal);
title('Gamma band');
xlabel('Time (s)');
ylabel('Voltage (uV)');

%% RÁCIO DE ENERGIAS ALPHA E GAMMA

% Set the interval duration in seconds
interval_duration = 2;

% Calculate the number of samples in each interval
interval_length = round(interval_duration * Fs);

% Calculate the number of intervals in the signal
num_intervals = floor(length(alpha_signal) / interval_length);

% Initialize an array to store the energy values
energy_alpha = zeros(num_intervals, 1);
energy_gamma = zeros(num_intervals, 1);

% Calculate the energy for each interval
for i = 1:num_intervals
    start_index = (i - 1) * interval_length + 1;
    end_index = start_index + interval_length - 1;
    energy_alpha(i) = sum(alpha_signal(start_index:end_index).^2);
    energy_gamma(i) = sum(gamma_signal(start_index:end_index).^2);
end

energy = energy_alpha ./ energy_gamma;

% Plot the energy values over time
t = (1:num_intervals) * interval_duration;
figure
plot(t, energy);
xlabel('Time (s)');
ylabel('Energy');
title('Signal Energy Over Time');

%% DASM - ALPHA

