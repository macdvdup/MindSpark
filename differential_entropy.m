function DE = differential_entropy(X, epoch_length, Fs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    DE = zeros(1,n_epochs-1);
    window = epoch_length * Fs;

    for i = 1:n_epochs
        epoch_signal = X((i-1)*window+1:i*window);
        sigma = var(epoch_signal);
        if i==1
            baseline = 0.5 * log(2 * pi * exp(1) * sigma);
        else
            DE(i-1) = 0.5 * log(2 * pi * exp(1) * sigma) - baseline;
        end
    end
end