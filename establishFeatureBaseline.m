function baseline = establishFeatureBaseline(signal)
    EalphaGamma_tp  = mean(signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./signal_energy(gamma_signaltp(m,:), epoch_length, Fs));
    EalphaBetaASM_tp  = mean(signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./asm_signal(vector_tp9,vector_tp10,epoch_length, Fs));
    EEbeta_tp  = mean(energy_entropy(beta_signaltp(m,:), epoch_length, Fs));
    EEgamma_tp  = mean(energy_entropy(gamma_signaltp(m,:), epoch_length, Fs));
    EEtheta_tp  = mean(energy_entropy(theta_signaltp(m,:), epoch_length,Fs));
    DASMalpha_tp  = mean(dasm_signal(vector_tp9(m,:), vector_tp10(m,:), epoch_length, Fs));
    PSD_alpha_tp  = mean(psd_signal(beta_signaltp(m,:), epoch_length, Fs, alpha_freq));
    
    EalphaGamma_af  = mean(signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./signal_energy(gamma_signalaf(m,:), epoch_length, Fs));
    EalphaBetaASM_af  = mean(signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./asm_signal(vector_af7,vector_af8,epoch_length, Fs));
    EEbeta_af  = mean(energy_entropy(beta_signalaf(m,:), epoch_length, Fs));
    EEgamma_af  = mean(energy_entropy(gamma_signalaf(m,:), epoch_length, Fs));
    EEtheta_af  = mean(energy_entropy(theta_signalaf(m,:), epoch_length,Fs));
    DASMalpha_af  = mean(dasm_signal(vector_af7(m,:), vector_af8(m,:), epoch_length, Fs));
    PSD_alpha_af  = mean(psd_signal(beta_signalaf(m,:), epoch_length, Fs, alpha_freq));
end