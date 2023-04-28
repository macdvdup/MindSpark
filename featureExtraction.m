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
alpha_bands = [8,12];
beta_bands = [12 30];
gamma_bands = [30 55];
bands = [alpha_bands; beta_bands; gamma_bands]; % alpha, beta, gamma

% Initialize matrices to store the feature values
E = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
%ER = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
EE = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
DE = zeros(n_epochs, size(bands,1), size(eeg_signal,2));
PSD = zeros(n_epochs, Fs+1,size(eeg_signal,2));
ASM = zeros(n_epochs, 2);
%RASM = zeros(n_epochs, size(eeg_signal,2));
DASM = zeros(n_epochs, 2);

% Loop over epochs and compute the features for each band and channel
for i = 1:n_epochs
    % Extract the current epoch
    epoch = eeg_signal((i-1)*epoch_samples+1 : i*epoch_samples,:);
    
    % Apply a window function to the epoch
    w = hamming(epoch_samples);
    epoch = epoch .* w;
    
    % Compute the power spectral density of the epoch
    [psd, f] = pwelch(epoch, w, epoch_samples/2, [], Fs);
    PSD(i,:) = psd;

    % Loop over the frequency bands and compute the features for each channel
    for j = 1:size(bands, 1)
        % Find the indices of the frequencies in the current band
        ind = (f >= bands(j,1)) & (f <= bands(j,2));
        % Compute the energy in the band for each channel
        energy = trapz(psd(ind,:));
        
        % Compute the total energy in the range of interest
        total_energy = trapz(psd(:,:), 1);
 
        % Compute the energy ratio for each channel
        %er = energy ./ total_energy; ER(i,j,:) = er;
        
        % Compute the energy entropy for each channel
        ee = -sum((energy ./ total_energy) .* log2(energy ./ total_energy), 1);
        
        % Store the feature values in the matrices
        E(i,j,:) = energy;
        EE(i,j,:) = ee;
        
        % Define parameters for differential entropy estimation
        nBins = 500; % number of bins for histogram
        for ch = k:size(epoch, 2)
            binEdges = linspace(min(epoch(:,k)), max(epoch(:,k)), nBins+1); % edges of histogram bins
            binWidth = binEdges(2) - binEdges(1); % width of histogram bins
            [binCounts, ~] = histcounts(eeg_signal(:,i), binEdges);
            binCounts = binCounts / sum(binCounts); % normalize bin counts to obtain PDF
            de = -sum(binCounts(binCounts > 0) .* log2(binCounts(binCounts > 0))) * binWidth;
            DE(i,j,k) = de;
        end  

        % Compute the asymmetry measures for the alpha band and each channel
        if j == 1 % alpha band
            left_channels = [1 2]; % TP9 and AF7
            right_channels = [3 4]; % TP10 and AF8
            
            % Compute the mean alpha energy in each hemisphere
            alpha_power_left = mean(energy(:,left_channels));
            alpha_power_right = mean(energy(:,right_channels));
            
            % Compute the asymmetry measures
            asm = alpha_power_right-alpha_power_left;
            % Store the asymmetry measures in the matrices
            ASM(i,:) = asm;
            DASM(i,:) = mean(DE(i,1,right_channels)) - mean(DE(i,1,left_channels))
        end
    end
end
