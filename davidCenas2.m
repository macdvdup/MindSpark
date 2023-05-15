clear; close all;
%addpath 'D:\Neuro\Projeto\MindSpark\functions'

% Load the EEG data and preprocess it as needed
eeg_signal_ini = csvread("my_matrix.csv");
Fs = 256; % Sampling frequency
eeg_signal_ini = eeg_signal_ini';

eeg_signal = eeg_signal_ini(:,1:4) - eeg_signal_ini(:,5);

% Pre-Processing channel signals
cutOffFreqs = [0.5, 80];
eeg_signal = preProcessEEG(eeg_signal, cutOffFreqs, Fs);

% Set the interval duration in seconds
interval_duration = 2;

% Calculate the number of samples in each interval
interval_length = round(interval_duration * Fs);


% Segment the EEG data into epochs
epoch_length = 2; % in seconds
epoch_samples = epoch_length * Fs; % fs is the sampling rate
n_epochs = floor(length(eeg_signal) / epoch_samples);

% Define the frequency bands of interest
bands = [8 12; 12 30; 30 80]; % alpha, beta, gamma

% Initialize matrices to store the feature values
E = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
ER = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
EE = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
DE = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
PSD = zeros(n_epochs, size(eeg_signal,2));
ASM = zeros(n_epochs, size(eeg_signal,2));
RASM = zeros(n_epochs, size(eeg_signal,2));
DASM = zeros(n_epochs, size(eeg_signal,2));

% Loop over epochs and compute the features for each band and channel
for i = 1:n_epochs
    % Extract the current epoch
    epoch = eeg_signal(:,(i-1)*epoch_samples+1 : i*epoch_samples);
    
    % Apply a window function to the epoch
    w = hamming(epoch_samples);
    epoch = epoch .* w;
    
    % Compute the power spectral density of the epoch
    [psd, f] = pwelch(epoch, w, epoch_samples/2, [], Fs);
    PSD(i,:,:) = psd;

    % Loop over the frequency bands and compute the features for each channel
    for j = 1:size(bands, 1)
        % Find the indices of the frequencies in the current band
        ind = (f >= bands(j,1)) & (f <= bands(j,2));
        % Compute the energy in the band for each channel
        energy = sum(psd(ind,:), 1);
        
        % Compute the total energy in the range of interest
        total_energy = sum(psd(:,:), 1);
        
        % Compute the energy ratio for each channel
        er = energy ./ total_energy;
        
        % Compute the energy entropy for each channel
        ee = -sum((energy ./ sum(energy)) .* log2(energy ./ sum(energy)), 1);
        
        % Compute the differential entropy for each channel
        de = -sum(   (psd(ind,:) ./ sum(psd(ind,:))      ) .* log2(psd(ind,:) ./ sum(psd(ind,:))), 1);
        
        % Store the feature values in the matrices
        E(i,j,:) = energy;
        ER(i,j,:) = er;
        EE(i,j,:) = ee;
        DE(i,j,:) = de;
        
        % Compute the asymmetry measures for the alpha band and each channel
        if j == 1 % alpha band
            left_channels = [1 2]; % TP9 and AF7
            right_channels = [3 4]; % TP10 and AF8
            
            % Compute the mean alpha power in each hemisphere
            alpha_power_left = mean(psd(ind,left_channels), 1);
            alpha_power_right = mean(psd(ind,right_channels), 1);
            
            % Compute the total alpha power for each hemisphere
            alpha_power_total_left = sum(psd(ind,left_channels), 1);
            alpha_power_total_right = sum(psd(ind,right_channels), 1);
            
            % Compute the asymmetry measures
            asm = (alpha_power_left - alpha_power_right) ./ (alpha_power_total_left + alpha_power_total_right);
            rasm = (alpha_power_left - alpha_power_right) ./ (alpha_power_left + alpha_power_right);
            dasm = abs(alpha_power_left - alpha_power_right);
            
            % Store the asymmetry measures in the matrices
            ASM(i,:) = asm;
            RASM(i,:) = rasm;
            DASM(i,:) = dasm;
        end
    end
end
