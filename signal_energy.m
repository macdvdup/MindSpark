% [Function] Compute the signal Energy
%
%     Developers: David Machado, Érica Gomes, José Almeida, Raquel Sousa
%     Neural Engineering Course, FEUP 2023
%
%     This function allows to obtain the signal energy for each epoch of
%     the input signal. This is one of the features mentioned as relevant
%     to evaluate the person state in the article:
%           Reference: Zhang, Y., Zhang, L., Hua, H., Jin, J., Zhu, L., 
%           Shu, L., Xu, X., Kuang, F., & Liu, Y. (2021). Relaxation Degree 
%           Analysis Using Frontal Electroencephalogram Under Virtual 
%           Reality Relaxation Scenes. Frontiers in Neuroscience. 
%                               https://doi.org/10.3389/fnins.2021.719869
%     
%      Input: X - Signal to process
%             epoch_length - The duration of each epoch
%             Fs - Sampling Frequency of the signal
%      Output: E - Vector with the energy of each epoch

function E = signal_energy(X, epoch_length, Fs)
    
    % Determines the number of points in each epoch
    epoch_samples = epoch_length * Fs;
    % Determines the number of epoch
    n_epochs = floor(length(X) / epoch_samples);
    % Creates the output vector
    E = zeros(1,n_epochs-1);
    % Stablishes the window based on Fs and duration of epoch
    window = epoch_length * Fs;
    
    % Computes the energy for each epoch
    for i = 1:n_epochs
        % Selects the epoch points
        epoch_signal = X((i-1)*window+1:i*window);
        % Determines the FFT of the epoch
        fx = fft(epoch_signal);
        
        % Determines the energy of each epoch using the methodology of the
        % above mentioned article
        if i==1
            % The first epoch serves as a baseline
            baseline = (1/(2*pi)) * sum(abs(fx).^2);
        else
            % Each epoch energy is then determined with reference to the
            % baseline
            E(i-1) = (1/(2*pi)) * sum(abs(fx).^2) - baseline;
        end
    end
end