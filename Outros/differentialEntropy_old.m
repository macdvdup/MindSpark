%% hjh
function [DE, da] = differentialEntropy(eeg_signal, epoch_length, sampling_freq)
    

    % Segment the EEG data into epochs
    epoch_samples = epoch_length * sampling_freq; % fs is the sampling rate
    n_epochs = floor(length(eeg_signal) / epoch_samples);
    
    % Define the frequency bands of interest
    delta_bands = [0 4];
    theta_bands = [4 8];
    alpha_bands = [8,13];
    beta_bands = [14 30];
    gamma_bands = [30 50];

    bands = [delta_bands; theta_bands; alpha_bands; beta_bands; gamma_bands];
    
    % Loop over epochs and compute the differential entropy for each band
    for i = 1:n_epochs
        % Extract the current epoch
        epoch = eeg_signal((i-1)*epoch_samples+1 : i*epoch_samples, :);
        
        % Apply a window function to the epoch
        w = hamming(epoch_samples);
        epoch = epoch .* w;
        
        % Compute the power spectral density of the epoch
        [psd, f] = pwelch(epoch, w, epoch_samples/2, [], fs);
        
        % Loop over the frequency bands and compute the differential entropy
        for j = 1:size(bands, 1)
            % Find the indices of the frequencies in the current band
            ind = (f >= bands(j,1)) & (f <= bands(j,2));
            
            % Compute the total power in the current band
            total_power = sum(psd(ind));
            
            % Normalize the power in the current band by the total power in the frequency range of interest
            relative_power = total_power / sum(psd(f>=0 & f<=80));
            
            % Compute the differential entropy for the current band
            DE(j,i) = -sum(relative_power .* log2(relative_power));
        end
    end
    
    % Compute the mean differential entropy across epochs for each band
    mean_H = mean(H, 2);
end

