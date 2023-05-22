function [alpha, theta, beta, gamma] = signal_decompose(signal, Fs, specific)
    alpha_freq = [8, 13];       
    theta_freq = [4, 8];
    beta_freq = [14, 30];
    gamma_freq = [30; 80];
    

    low_cf = 1;
    high_cf = 80;
    filter_order = 4;
    [b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
    signal = filtfilt(b,a,signal);
    
    d = designfilt('bandstopiir','FilterOrder',4, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',Fs);
    signal = filtfilt(d,signal);

    if ~exist('specific', 'var')
        alpha = eegfilt(signal, Fs, alpha_freq(1), alpha_freq(2));
        beta = eegfilt(signal, Fs, beta_freq(1), beta_freq(2));
        theta = eegfilt(signal, Fs, theta_freq(1), theta_freq(2));
        gamma = eegfilt(signal, Fs, gamma_freq(1), gamma_freq(2));
    elseif strcmp(specific, 'beta')
        beta = eegfilt(signal, Fs, beta_freq(1), beta_freq(2));
        alpha = []; theta = []; gamma = [];
    elseif strcmp(specific, 'alpha')
        alpha = eegfilt(signal, Fs, alpha_freq(1), alpha_freq(2));
        beta = []; theta = []; gamma = [];
    end
end