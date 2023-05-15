function E = signal_energy(X, epoch_length, Fs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    E = zeros(1,n_epochs);
    window = epoch_length * Fs;

    for i = 1:n_epochs
        epoch_signal = X((i-1)*window+1:i*window);
        fx = fft(epoch_signal);
        E(i) = (1/(2*pi)) * sum(abs(fx).^2);
    end
end