clc; clear all; close all;

Fs = 128;

% Get the mean signal for relax state
subjects = dir('GAMEEMO');

% CALM
calm_AF3 = [];
for i=1:3
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF3 = [calm_AF3; AF3'];
    end
end
calm_AF3 = mean(calm_AF3, 1);

calm_AF4 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF4 = [calm_AF4; AF4'];
    end
end
calm_AF4 = mean(calm_AF4, 1);

% HORROR
horror_AF3 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF3 = [horror_AF3; AF3'];
    end
end
horror_AF3 = mean(horror_AF3, 1);

horror_AF4 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF4 = [horror_AF4; AF4'];
    end
end
horror_AF4 = mean(horror_AF4, 1);

low_cf = 1;
high_cf = 50;
filter_order = 4;
[b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));

calm_AF3 = filtfilt(b,a,calm_AF3);
horror_AF3 = filtfilt(b,a,horror_AF3);
calm_AF4 = filtfilt(b,a,calm_AF4);
horror_AF4 = filtfilt(b,a,horror_AF4);
% 
% Wo = 50/(Fs/2);  BW = Wo/35;
% [bn,an] = iirnotch(Wo,BW);
% 
% calm_AF3 = filtfilt(bn,an,calm_AF3);
% horror_AF3 = filtfilt(bn,an,horror_AF3);
% calm_AF4 = filtfilt(bn,an,calm_AF4);
% horror_AF4 = filtfilt(bn,an,horror_AF4);

eeg_signal_AF3 = [calm_AF3; horror_AF3];
eeg_signal_AF4 = [calm_AF4; horror_AF4];

% Define frequency bands of interest
alpha_freq = [8, 13];
theta_freq = [4, 8];
gamma_freq = [30; 50];

% Apply bandpass filter to signal
alpha_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, alpha_freq(1), alpha_freq(2));
theta_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, theta_freq(1), theta_freq(2));
gamma_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, gamma_freq(1), gamma_freq(2));
alpha_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, alpha_freq(1), alpha_freq(2));
theta_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, theta_freq(1), theta_freq(2));
gamma_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, gamma_freq(1), gamma_freq(2));

% ALPHA POWER SPECTRAL DENSITY (PSD)

window = 30*Fs;
alpha_power = zeros(size(eeg_signal_AF3,1), length(eeg_signal_AF3)-window+1);

for k = 1:size(alpha_signal_AF3,1)
    channel_signal_AF3 = alpha_signal_AF3(k,:);
    channel_signal_AF4 = alpha_signal_AF4(k,:);
    for i = 1:length(alpha_power)
        signal_window_AF3 = channel_signal_AF3(i:i+window-1);
        signal_window_AF4 = channel_signal_AF4(i:i+window-1);
    
        % Compute power spectral density
        [psd_AF3, freq_AF3] = pwelch(signal_window_AF3, [], [], [], Fs);
        [psd_AF4, freq_AF4] = pwelch(signal_window_AF4, [], [], [], Fs);

        % Integrate PSD over alpha frequency range to obtain alpha power
        alpha_freq_indices_AF3 = freq_AF3 >= alpha_freq(1) & freq_AF3 <= alpha_freq(2);
        alpha_freq_indices_AF4 = freq_AF4 >= alpha_freq(1) & freq_AF4 <= alpha_freq(2);

        alpha_power(k, i) = (trapz(freq_AF3(alpha_freq_indices_AF3), psd_AF3(alpha_freq_indices_AF3))+trapz(freq_AF4(alpha_freq_indices_AF4), psd_AF4(alpha_freq_indices_AF4)))./2;
    end
end

figure
t = linspace(0, length(eeg_signal_AF3)/Fs, length(alpha_power(1,:)));
plot(t, alpha_power')
title('Alpha Power')
drawnow


% THETA POWER SPECTRAL DENSITY (PSD)

window = 30*Fs;
theta_power = zeros(size(eeg_signal_AF3,1), length(eeg_signal_AF3)-window+1);

for k = 1:size(theta_signal_AF3,1)
    channel_signal_AF3 = theta_signal_AF3(k,:);
    channel_signal_AF4 = theta_signal_AF4(k,:);
    for i = 1:length(theta_power)
        signal_window_AF3 = channel_signal_AF3(i:i+window-1);
        signal_window_AF4 = channel_signal_AF4(i:i+window-1);
    
        % Compute power spectral density
        [psd_AF3, freq_AF3] = pwelch(signal_window_AF3, [], [], [], Fs);
        [psd_AF4, freq_AF4] = pwelch(signal_window_AF4, [], [], [], Fs);

        % Integrate PSD over alpha frequency range to obtain alpha power
        theta_freq_indices_AF3 = freq_AF3 >= theta_freq(1) & freq_AF3 <= theta_freq(2);
        theta_freq_indices_AF4 = freq_AF4 >= theta_freq(1) & freq_AF4 <= theta_freq(2);

        theta_power(k, i) = (trapz(freq_AF3(theta_freq_indices_AF3), psd_AF3(theta_freq_indices_AF3))+trapz(freq_AF4(theta_freq_indices_AF4), psd_AF4(theta_freq_indices_AF4)))./2;
    end
end
figure
t = linspace(0, length(eeg_signal_AF3)/Fs,  length(theta_power(1,:)));
plot(t, theta_power')
title('Theta Power')
drawnow

% GAMMA POWER SPECTRAL DENSITY (PSD)

window = 30*Fs;
gamma_power = zeros(size(eeg_signal_AF3,1), length(eeg_signal_AF3)-window+1);

for k = 1:size(gamma_signal_AF3,1)
    channel_signal_AF3 = gamma_signal_AF3(k,:);
    channel_signal_AF4 = gamma_signal_AF4(k,:);
    for i = 1:length(gamma_power)
        signal_window_AF3 = channel_signal_AF3(i:i+window-1);
        signal_window_AF4 = channel_signal_AF4(i:i+window-1);
    
        % Compute power spectral density
        [psd_AF3, freq_AF3] = pwelch(signal_window_AF3, [], [], [], Fs);
        [psd_AF4, freq_AF4] = pwelch(signal_window_AF4, [], [], [], Fs);

        % Integrate PSD over alpha frequency range to obtain alpha power
        gamma_freq_indices_AF3 = freq_AF3 >= gamma_freq(1) & freq_AF3 <= gamma_freq(2);
        gamma_freq_indices_AF4 = freq_AF4 >= gamma_freq(1) & freq_AF4 <= gamma_freq(2);

        gamma_power(k, i) = (trapz(freq_AF3(gamma_freq_indices_AF3), psd_AF3(gamma_freq_indices_AF3))+trapz(freq_AF4(gamma_freq_indices_AF4), psd_AF4(gamma_freq_indices_AF4)))./2;
    end
end
figure
t = linspace(0, length(eeg_signal_AF3)/Fs,  length(gamma_power(1,:)));
plot(t, gamma_power')
title('Gamma Power')
drawnow

% RELATIVE POWER SPECTRAL DENSITY (PSD)

window = 30*Fs;
relative_power = alpha_power./(alpha_power+theta_power+gamma_power);

figure
t = linspace(0, length(eeg_signal_AF3)/Fs,  length(relative_power(1,:)));
plot(t, relative_power')
title('Relative Alpha Power')
drawnow

%% TAKING TOTAL POWER INTO ACCOUNT




%% TEST DASM
clc; clear all;
fs = 128;

% Get the mean signal for relax state
subjects = dir('GAMEEMO');

% CALM
calm_AF3 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF3 = [calm_AF3; AF3'];
    end
end
calm_AF3 = mean(calm_AF3, 1);

calm_AF4 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF4 = [calm_AF4; AF4'];
    end
end
calm_AF4 = mean(calm_AF4, 1);

% HORROR
horror_AF3 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF3 = [horror_AF3; AF3'];
    end
end
horror_AF3 = mean(horror_AF3, 1);

horror_AF4 = [];
for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF4 = [horror_AF4; AF4'];
    end
end
horror_AF4 = mean(horror_AF4, 1);

low_cf = 1;
high_cf = 35;
filter_order = 6;
[b, a] = butter(filter_order, [low_cf high_cf]/(fs/2));

calm_AF3 = filtfilt(b,a,calm_AF3);
calm_AF4 = filtfilt(b,a,calm_AF4);
horror_AF3 = filtfilt(b,a,horror_AF3);
horror_AF4 = filtfilt(b,a,horror_AF4);

f_alpha = [8 13];

epoch_length = 30;
epoch_samples = epoch_length * fs; % fs is the sampling rate
n_epochs = floor(length(calm_AF3) / epoch_samples);

AF3 = [calm_AF3; horror_AF3];
AF4 = [calm_AF4; horror_AF4];

differential_entropy = [];
for k = 1:size(AF3,1)
    de_vector = [];
    for i = 1:n_epochs
        % Extract the current epoch
        epoch_AF3 = AF3(k,(i-1)*epoch_samples+1 : i*epoch_samples)';
        epoch_AF4 = AF4(k,(i-1)*epoch_samples+1 : i*epoch_samples)';
        
        % Apply a window function to the epoch
        w = hamming(epoch_samples);
        epoch_AF3 = epoch_AF3 .* w;
        epoch_AF4 = epoch_AF4 .* w;
    
        % Compute the power spectral density of the epoch
        [psd_tp9, freq_tp9] = pwelch(epoch_AF3, w, epoch_samples/2, [], fs);
        [psd_tp10, freq_tp10] = pwelch(epoch_AF4, w, epoch_samples/2, [], fs);
        
        % Calculate the spectral entropy of the PSD
        alpha_indices_tp9 = freq_tp10 >= f_alpha(1) & freq_tp9 <= f_alpha(2); % indices of alpha band
        alpha_indices_tp10 = freq_tp10 >= f_alpha(1) & freq_tp9 <= f_alpha(2); % indices of alpha band
        psd_alpha_tp9 = mean(psd_tp9(alpha_indices_tp9, :, :), 1); % average PSD in alpha band
        psd_alpha_tp10 = mean(psd_tp10(alpha_indices_tp10, :, :), 1);
        spectral_entropy_tp9 = -sum(psd_alpha_tp9 .* log2(psd_alpha_tp9), 1); % spectral entropy
        spectral_entropy_tp10 = -sum(psd_alpha_tp10 .* log2(psd_alpha_tp10), 1);
    
        % Convert spectral entropy to differential entropy
        value_tp9 = spectral_entropy_tp9 ./ log2(exp(1));
        value_tp10 = spectral_entropy_tp10 ./ log2(exp(1));
        
        de_vector(end+1) = value_tp9 - value_tp10;
    end
    differential_entropy = [differential_entropy; de_vector];
end

% Plot the differential entropy
figure;
plot(differential_entropy');
xlabel('Segment number');
ylabel('Differential entropy');
title('DAMS para o AF3 e AF4');

