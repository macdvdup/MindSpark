clc; clear all; close all;
addpath 'D:\Neuro\Projeto\MindSpark\functions\sigprocfunc'

mu = MuseUdp();

epoch_length = 30;
Fs = 256;

tic
cont = true;
vector_tp9 = [];
vector_af7 = [];
vector_af8 = [];
vector_tp10 = [];

EalphaGamma_tp = [];
EalphaBetaASM_tp = [];
EEbeta_tp = [];
EEgamma_tp = [];
EEtheta_tp = [];
DASMalpha_tp = [];
PSD_alpha_tp = [];

EalphaGamma_af = [];
EalphaBetaASM_af = [];
EEbeta_af = [];
EEgamma_af = [];
EEtheta_af = [];
DASMalpha_af = [];
PSD_alpha_af = [];

alpha_freq = [8, 13];       
theta_freq = [4, 8];
beta_freq = [14, 30];
gamma_freq = [30; 50];


while cont
    [data, timestamp, success] = mu.get_eeg_sample();
    if success
        vector_tp9(end+1) = data(1)-data(5);
        vector_af7(end+1) = data(2)-data(5);
        vector_af8(end+1) = data(3)-data(5);
        vector_tp10(end+1) = data(4)-data(5);
    end

    if toc > 30
        cont = false;
    end
end

% PROCESSAMENTO E EXTRAÇÃO DE FEATURES
        
signal_tp = (vector_tp9 + vector_tp10) ./2;
signal_af = (vector_af8 + vector_af7) ./2;

[alpha_signaltp, theta_signaltp, beta_signaltp, gamma_signaltp] = signal_decompose(signal_tp, Fs);
[alpha_signalaf, theta_signalaf, beta_signalaf, gamma_signalaf] = signal_decompose(signal_af, Fs);

for m=1:size(alpha_signalaf)
    EalphaGamma_tp(end+1) = mean(signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./signal_energy(gamma_signaltp(m,:), epoch_length, Fs));
    EalphaBetaASM_tp(end+1) = mean(signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./asm_signal(vector_tp9,vector_tp10,epoch_length, Fs));
    EEbeta_tp(end+1) = mean(energy_entropy(beta_signaltp(m,:), epoch_length, Fs));
    EEgamma_tp(end+1) = mean(energy_entropy(gamma_signaltp(m,:), epoch_length, Fs));
    EEtheta_tp(end+1) = mean(energy_entropy(theta_signaltp(m,:), epoch_length,Fs));
    DASMalpha_tp(end+1) = mean(dasm_signal(vector_tp9(m,:), vector_tp10(m,:), epoch_length, Fs));
    PSD_alpha_tp(end+1) = mean(psd_signal(beta_signaltp(m,:), epoch_length, Fs, alpha_freq));
    
    EalphaGamma_af(end+1) = mean(signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./signal_energy(gamma_signalaf(m,:), epoch_length, Fs));
    EalphaBetaASM_af(end+1) = mean(signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./asm_signal(vector_af7,vector_af8,epoch_length, Fs));
    EEbeta_af(end+1) = mean(energy_entropy(beta_signalaf(m,:), epoch_length, Fs));
    EEgamma_af(end+1) = mean(energy_entropy(gamma_signalaf(m,:), epoch_length, Fs));
    EEtheta_af(end+1) = mean(energy_entropy(theta_signalaf(m,:), epoch_length,Fs));
    DASMalpha_af(end+1) = mean(dasm_signal(vector_af7(m,:), vector_af8(m,:), epoch_length, Fs));
    PSD_alpha_af(end+1) = mean(psd_signal(beta_signalaf(m,:), epoch_length, Fs, alpha_freq));
end