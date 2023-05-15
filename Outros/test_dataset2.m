clc; clear all; close all;

Fs = 128;
addpath 'functions\'
% Get the mean signal for relax state
subjects = dir('GAMEEMO');

patients_results = zeros(2, 2);

for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF3 = AF3';
        calm_AF4 = AF4';

        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF3 = AF3';
        horror_AF4 = AF4';

        low_cf = 1;
        high_cf = 50;
        filter_order = 4;
        [b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));

        calm_AF3 = filtfilt(b,a,calm_AF3);
        horror_AF3 = filtfilt(b,a,horror_AF3);
        
        eeg_signal_AF3 = [calm_AF3; horror_AF3];
        eeg_signal_AF4 = [calm_AF4; horror_AF4];

        alpha_freq = [8, 13];       
        theta_freq = [4, 8];
        beta_freq = [14, 30];
        gamma_freq = [30; 50];

        alpha_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, beta_freq(1), beta_freq(2));
        theta_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, theta_freq(1), theta_freq(2));
        gamma_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, gamma_freq(1), gamma_freq(2));
        alpha_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, beta_freq(1), beta_freq(2));
        theta_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, theta_freq(1), theta_freq(2));
        gamma_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, gamma_freq(1), gamma_freq(2));

        % GAMMA POWER SPECTRAL DENSITY (PSD)
        
        epoch_length = 30; % in seconds
        epoch_samples = epoch_length * Fs; % fs is the sampling rate
        n_epochs = floor(length(eeg_signal_AF3) / epoch_samples);

        window = epoch_length*Fs;
        gamma_power = [];
        alpha_power = [];
        beta_power = [];
        theta_power = [];

        for k = 1:size(gamma_signal_AF3,1)
            channel_signal_AF3g = gamma_signal_AF3(k,:);
            channel_signal_AF4g = gamma_signal_AF4(k,:);

            channel_signal_AF3a = alpha_signal_AF3(k,:);
            channel_signal_AF4a = alpha_signal_AF4(k,:);

            channel_signal_AF3b = beta_signal_AF3(k,:);
            channel_signal_AF4b = beta_signal_AF4(k,:);

            channel_signal_AF3t = theta_signal_AF3(k,:);
            channel_signal_AF4t = theta_signal_AF4(k,:);
            for m = 1:n_epochs
                signal_window_AF3g = channel_signal_AF3g(m:m+window-1);
                signal_window_AF4g = channel_signal_AF4g(m:m+window-1);

                signal_window_AF3a = channel_signal_AF3a(m:m+window-1);
                signal_window_AF4a = channel_signal_AF4a(m:m+window-1);
                
                signal_window_AF3b = channel_signal_AF3b(m:m+window-1);
                signal_window_AF4b = channel_signal_AF4b(m:m+window-1);

                signal_window_AF3t = channel_signal_AF3t(m:m+window-1);
                signal_window_AF4t = channel_signal_AF4t(m:m+window-1);

%                w = hamming(window);
%                signal_window_AF3 = signal_window_AF3 .* w;
%                signal_window_AF4 = signal_window_AF4 .* w;

                % Compute power spectral density
                [psd_AF3g, freq_AF3g] = pwelch(signal_window_AF3g, [], [], [], Fs);
                [psd_AF4g, freq_AF4g] = pwelch(signal_window_AF4g, [], [], [], Fs);
                
                [psd_AF3a, freq_AF3a] = pwelch(signal_window_AF3a, [], [], [], Fs);
                [psd_AF4a, freq_AF4a] = pwelch(signal_window_AF4a, [], [], [], Fs);

                [psd_AF3b, freq_AF3b] = pwelch(signal_window_AF3g, [], [], [], Fs);
                [psd_AF4b, freq_AF4b] = pwelch(signal_window_AF4g, [], [], [], Fs);
                
                [psd_AF3t, freq_AF3t] = pwelch(signal_window_AF3t, [], [], [], Fs);
                [psd_AF4t, freq_AF4t] = pwelch(signal_window_AF4t, [], [], [], Fs);


                % Integrate PSD over alpha frequency range to obtain alpha power
                gamma_freq_indices_AF3 = freq_AF3g >= gamma_freq(1) & freq_AF3g <= gamma_freq(2);
                gamma_freq_indices_AF4 = freq_AF4g >= gamma_freq(1) & freq_AF4g <= gamma_freq(2);
                
                alpha_freq_indices_AF3 = freq_AF3a >= alpha_freq(1) & freq_AF3a <= alpha_freq(2);
                alpha_freq_indices_AF4 = freq_AF4a >= alpha_freq(1) & freq_AF4a <= alpha_freq(2);
                
                beta_freq_indices_AF3 = freq_AF3b >= beta_freq(1) & freq_AF3b <= beta_freq(2);
                beta_freq_indices_AF4 = freq_AF4b >= beta_freq(1) & freq_AF4b <= beta_freq(2);
                
                theta_freq_indices_AF3 = freq_AF3t >= theta_freq(1) & freq_AF3t <= theta_freq(2);
                theta_freq_indices_AF4 = freq_AF4t >= theta_freq(1) & freq_AF4t <= theta_freq(2);

                gamma_power(k, m) = (trapz(freq_AF3g(gamma_freq_indices_AF3), psd_AF3g(gamma_freq_indices_AF3))+trapz(freq_AF4g(gamma_freq_indices_AF4), psd_AF4g(gamma_freq_indices_AF4)))./2;
                alpha_power(k, m) = (trapz(freq_AF3a(alpha_freq_indices_AF3), psd_AF3a(alpha_freq_indices_AF3))+trapz(freq_AF4a(alpha_freq_indices_AF4), psd_AF4a(alpha_freq_indices_AF4)))./2;
                beta_power(k, m) = (trapz(freq_AF3b(beta_freq_indices_AF3), psd_AF3b(beta_freq_indices_AF3))+trapz(freq_AF4b(beta_freq_indices_AF4), psd_AF4b(beta_freq_indices_AF4)))./2;
                theta_power(k, m) = (trapz(freq_AF3t(theta_freq_indices_AF3), psd_AF3t(theta_freq_indices_AF3))+trapz(freq_AF4t(theta_freq_indices_AF4), psd_AF4t(theta_freq_indices_AF4)))./2;
            end
        end
        %patients_results(1:2, i) = mean(beta_power./(gamma_power + theta_power + beta_power) , 2);
        patients_results(1:2, i) = mean(alpha_power, 2);
        %patients_results(i, :) = [];
    end
end

t = linspace(1,length(patients_results), length(patients_results));
plot(t, patients_results);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('PSD');
title('PSD Theta Signal');


