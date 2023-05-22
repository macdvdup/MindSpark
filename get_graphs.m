addpath 'functions\'
% Get the mean signal for relax state
% subjects = dir('GAMEEMO');
% 
% sub = 7;
% 
% file_str = strcat(replace(subjects(sub).name, {'(', ')'}, ''), 'G2AllRawChannels.mat');
% str = {'GAMEEMO', subjects(sub).name, 'Raw EEG Data', '.mat format', file_str};
% str = strjoin(str, '\');
% load(str);
% eeg_data = AF3';
% Fs = 128;
thenvelopewindow=20;
eeg_data = csvread("Outros/my_matrix.csv");
Fs = 256; % Sampling frequency
eeg_data = eeg_data(1,10*Fs:Fs*21);
t = linspace(0, length(eeg_data)/Fs,length(eeg_data));

T = 1/Fs; % Sampling period
L = length(eeg_data); % Length of signal
Y = fft(eeg_data); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
[P1,maxChunkTime] = envelope(P1, thenvelopewindow, 'rms');
f = Fs*(0:(L/2))/L; % Frequency vector
figure
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Raw Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 80])

figure;
plot(t, eeg_data);
title('Raw EEG signal');
xlabel('Time (s)');
ylabel('Voltage (uV)');

low_cf = 1;
high_cf = 80;
filter_order = 4;
[b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
filtered_signal = filtfilt(b,a,eeg_data);

T = 1/Fs; % Sampling period
L = length(filtered_signal); % Length of signal
Y = fft(filtered_signal); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L; % Frequency vector
figure
[P1,maxChunkTime] = envelope(P1, thenvelopewindow, 'rms');
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Filtered Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 80])
figure;
plot(t,filtered_signal);
title('Filtered EEG signal');
xlabel('Time (s)');
ylabel('Voltage (uV)');

d = designfilt('bandstopiir','FilterOrder',4, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',Fs);
notch_signal = filtfilt(d,filtered_signal);

figure;
plot(t,notch_signal);
title('EEG Signal with Notch Filter');
xlabel('Time (s)');
ylabel('Voltage (uV)');

T = 1/Fs; % Sampling period
L = length(notch_signal); % Length of signal
Y = fft(notch_signal); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
[P1,maxChunkTime] = envelope(P1, thenvelopewindow, 'rms');
f = Fs*(0:(L/2))/L; % Frequency vector
figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Filtered-Notch Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 80])
% Define frequency bands of interest
alpha_freq = [8, 13];
theta_freq = [4, 8];
beta_freq = [13, 30];
gamma_freq = [30, 100];

[alpha_signal, theta_signal, beta_signal, gamma_signal] = signal_decompose(notch_signal, Fs);

figure;
subplot(4,1,1);
plot(t,theta_signal);
title('Theta band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,2);
plot(t,alpha_signal);
title('Alpha band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,3);
plot(t,beta_signal);
title('Beta band');
xlabel('Time (s)');
ylabel('Voltage (uV)');
subplot(4,1,4);
plot(t,gamma_signal);
title('Gamma band');
xlabel('Time (s)');
ylabel('Voltage (uV)');