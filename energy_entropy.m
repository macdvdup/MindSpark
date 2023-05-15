function EE = energy_entropy(X, epoch_length, Fs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    window = epoch_length * Fs;
    EE = zeros(1,n_epochs);

    for i = 1:n_epochs
        epoch_signal = X((i-1)*window+1:i*window);
        [psd, ~] = pwelch(epoch_signal, [], [], [], Fs);
        psd = psd / sum(psd);
        EE(i) = -sum(log(psd.^2));
    end
end