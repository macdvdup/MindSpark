function DE = differential_entropy(X, epoch_length, Fs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    DE = zeros(1,n_epochs);
    window = epoch_length * Fs;

    for i = 1:n_epochs
        epoch_signal = X(i:i+window-1);
        sigma = var(epoch_signal);
        DE(i) = 0.5 * log(2 * pi * exp(1) * sigma);
    end
end