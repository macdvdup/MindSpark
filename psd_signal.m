function PSD = psd_signal(X, epoch_length, Fs, Freqs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    PSD = zeros(1,n_epochs);
    window = epoch_length * Fs;

    for i = 1:n_epochs
        epoch_signal = X(i:i+window-1);
        [psd, f] = pwelch(epoch_signal, [], [], [], Fs);
        ind = f >= Freqs(1) & f <= Freqs(2);
        PSD(i) = (trapz(f(ind), psd(ind)));
    end
end