clear all; close all;

addpath 'functions'

data = csvread("my_matrix.csv");
fs = 256; % Sampling frequency
data = data(:,1:4096);

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
[b, a] = butter(filter_order, [low_cf high_cf]/(fs/2));

vector_tp9 = filtfilt(b,a,vector_tp9);
vector_tp9 = vector_tp9';

vector_tp10 = filtfilt(b,a,vector_tp10);
vector_tp10 = vector_tp10';

f_alpha = [8 13];

epoch_length = 2;
epoch_samples = epoch_length * fs; % fs is the sampling rate
n_epochs = floor(length(vector_tp9) / epoch_samples);

differential_entropy = [];
for i = 1:n_epochs
    % Extract the current epoch
    epoch_tp9 = vector_tp9((i-1)*epoch_samples+1 : i*epoch_samples, :);
    epoch_tp10 = vector_tp10((i-1)*epoch_samples+1 : i*epoch_samples, :);

    % Apply a window function to the epoch
    w = hamming(epoch_samples);
    epoch_tp9 = epoch_tp9 .* w;
    epoch_tp10 = epoch_tp10 .* w;

    % Compute the power spectral density of the epoch
    [psd_tp9, freq_tp9] = pwelch(epoch_tp9, w, epoch_samples/2, [], fs);
    [psd_tp10, freq_tp10] = pwelch(epoch_tp10, w, epoch_samples/2, [], fs);
    
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
    
    differential_entropy(end+1) = value_tp9 - value_tp10;
end


% Plot the differential entropy
figure;
plot(differential_entropy);
xlabel('Segment number');
ylabel('Differential entropy');
title('DAMS para o TP9 e TP10');
