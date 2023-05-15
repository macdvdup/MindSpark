function DE = differentialEntropy(signal, epoch_length, Fs)

    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(eeg_signal_AF3) / epoch_samples);
    
    for m = 1:n_epochs
        signal_window = signal;
        
    end

end