clc; clear all; close all;
addpath 'D:\Neuro\Projeto\MindSpark\functions\sigprocfunc'

mu = MuseUdp();

epoch_length = 30;
Fs = 256;

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

% Establish baseline
while cont
    [data, timestamp, success] = mu.get_eeg_sample();
    if success
        vector_tp9(end+1) = data(1)-data(5);
        vector_af7(end+1) = data(2)-data(5);
        vector_af8(end+1) = data(3)-data(5);
        vector_tp10(end+1) = data(4)-data(5);
    end

    if flagSessionEnded
        cont = false;
    end
end

% PROCESSAMENTO E EXTRAÇÃO DE FEATURES
        
signal_tp = (vector_tp9 + vector_tp10) ./2;
signal_af = (vector_af8 + vector_af7) ./2;

[alpha_signaltp, theta_signaltp, beta_signaltp, gamma_signaltp] = signal_decompose(signal_tp, Fs);
[alpha_signalaf, theta_signalaf, beta_signalaf, gamma_signalaf] = signal_decompose(signal_af, Fs);

EalphaGamma_af = (signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./signal_energy(gamma_signalaf(m,:), epoch_length, Fs));
EalphaBetaASM_af = (signal_energy(alpha_signalaf(m,:), epoch_length, Fs)./asm_signal(vector_af7,vector_af8,epoch_length, Fs));
EEbeta_af = (energy_entropy(beta_signalaf(m,:), epoch_length, Fs));
EEgamma_af = (energy_entropy(gamma_signalaf(m,:), epoch_length, Fs));
EEtheta_af = (energy_entropy(theta_signalaf(m,:), epoch_length,Fs));
DASMalpha_af = (dasm_signal(vector_af7(m,:), vector_af8(m,:), epoch_length, Fs));
PSD_beta_af = (psd_signal(beta_signalaf(m,:), epoch_length, Fs, alpha_freq));

EalphaGamma_tp = (signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./signal_energy(gamma_signaltp(m,:), epoch_length, Fs));
EalphaBetaASM_tp = (signal_energy(alpha_signaltp(m,:), epoch_length, Fs)./asm_signal(vector_tp9,vector_tp10,epoch_length, Fs));
EEbeta_tp = (energy_entropy(beta_signaltp(m,:), epoch_length, Fs));
EEgamma_tp = (energy_entropy(gamma_signaltp(m,:), epoch_length, Fs));
EEtheta_tp = (energy_entropy(theta_signaltp(m,:), epoch_length,Fs));
DASMalpha_tp = (dasm_signal(vector_tp9(m,:), vector_tp10(m,:), epoch_length, Fs));
PSD_beta_tp = (psd_signal(beta_signaltp(m,:), epoch_length, Fs, alpha_freq));


datasetFolder = {pwd, "creating dataset"};
modelFile = strjoin(datasetFolder, '\')+"\bestModelBoth.mat";
load(modelFile,'bestClassifier','mu','sigma');

features = [EalphaGamma_af' EalphaBetaASM_af' EEbeta_af' EEgamma_af' EEtheta_af' DASMalpha_af' PSD_beta_af' EalphaGamma_tp' EalphaBetaASM_tp' EEbeta_tp' EEgamma_tp' EEtheta_tp' DASMalpha_tp' PSD_beta_tp'];
features = (features-mu)./sigma;
[predictedY, scores, as] = predict(bestClassifier, X_test);
% predictedY- 1 if stressed, 0 if normal
% scores - analog to class probability

for i=1:length(predictedY)
    str=""
    if predictedY(i)==1
        str = str + "Stress: Level " + int2str(scores(i)*5)
    else
        str = str + "Calm : Level " + int2str(scores(i)*5)
    end
end

plot(scores)
