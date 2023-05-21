function PSD = psd_signal(X, epoch_length, Fs, Freqs)
    epoch_samples = epoch_length * Fs;
    n_epochs = floor(length(X) / epoch_samples);
    PSD = zeros(1,n_epochs-1);
    window = epoch_length * Fs;
    
    for i = 1:n_epochs
        epoch_signal = X((i-1)*window+1:i*window);
        [psd, f] = pwelch(epoch_signal, [], [], [], Fs);
        ind = f >= Freqs(1) & f <= Freqs(2);
        if i==1  
            baseline = (trapz(f(ind), psd(ind))/trapz(f, psd));
        else
            PSD(i-1) = (trapz(f(ind), psd(ind))/trapz(f, psd))-baseline;
        end
    end
end