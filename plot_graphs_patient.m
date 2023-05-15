clc; clear all; close all;

Fs = 128;
addpath 'D:\Neuro\Projeto\MindSpark\functions\sigprocfunc'
% Get the mean signal for relax state
subjects = dir('GAMEEMO');

epoch_length = 30; % in seconds
psd_patients = [];
e_patients = [];
ee_patients = [];
dasm_patients = [];
de_patients = [];

EalphaGamma = [];
EalphaBetaASM = [];
EEbeta = [];
EEgamma = [];
EEtheta = [];
DASMalpha = [];

for i=2:2
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF3 = AF3';
        calm_AF4 = AF4';

        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF3 = AF3';
        horror_AF4 = AF4';

        low_cf = 1;
        high_cf = 50;
        filter_order = 4;
        [b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));

        calm_AF3 = filtfilt(b,a,calm_AF3);
        horror_AF3 = filtfilt(b,a,horror_AF3);
        
        eeg_signal_AF3 = [calm_AF3; horror_AF3];
        eeg_signal_AF4 = [calm_AF4; horror_AF4];

        alpha_freq = [8, 13];       
        theta_freq = [4, 8];
        beta_freq = [14, 30];
        gamma_freq = [30; 50];

        alpha_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, beta_freq(1), beta_freq(2));
        theta_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, theta_freq(1), theta_freq(2));
        gamma_signal_AF3 = eegfilt(eeg_signal_AF3, Fs, gamma_freq(1), gamma_freq(2));
        alpha_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, beta_freq(1), beta_freq(2));
        theta_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, theta_freq(1), theta_freq(2));
        gamma_signal_AF4 = eegfilt(eeg_signal_AF4, Fs, gamma_freq(1), gamma_freq(2));
        
        alpha_signal = (alpha_signal_AF4 + alpha_signal_AF3) ./2;
        beta_signal = (beta_signal_AF4 + beta_signal_AF3) ./2;
        gamma_signal = (gamma_signal_AF4 + gamma_signal_AF3) ./2;
        theta_signal = (theta_signal_AF4 + theta_signal_AF3) ./2;

        for m=1:size(eeg_signal_AF3)
            EalphaGamma(m,:) = signal_energy(alpha_signal(m,:), epoch_length, Fs)./signal_energy(gamma_signal(m,:), epoch_length, Fs);
            EalphaBetaASM(m,:) = signal_energy(alpha_signal(m,:), epoch_length, Fs)./asm_signal(beta_signal_AF3,beta_signal_AF4,epoch_length, Fs);
            EEbeta(m,:) = energy_entropy(beta_signal(m,:), epoch_length, Fs);
            EEgamma(m,:) = energy_entropy(gamma_signal(m,:), epoch_length, Fs);
            EEtheta(m,:) = energy_entropy(theta_signal(m,:), epoch_length,Fs);
            DASM(m,:) = dasm_signal(alpha_signal_AF3(m,:), alpha_signal_AF4(m,:), epoch_length, Fs);
            PSD(m,:) = psd_signal(beta_signal(m,:), epoch_length, Fs, alpha_freq);
        end
    end
end

t = linspace(1,length(calm_AF3)/Fs, length(EalphaGamma));
figure
plot(t, EalphaGamma);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('E alpha/E gamma - P7/P8');

figure
plot(t, EalphaBetaASM);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('E alpha/ASM beta - P7/P8');

figure
plot(t, EEbeta);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('EE beta - P7/P8');

figure
plot(t, EEgamma);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('EE gamma - P7/P8');

figure
plot(t, EEtheta);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('EE theta - P7/P8');

figure
plot(t, DASM);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('DASM alpha - P7/P8');

figure
plot(t, PSD);
legend("Calm","Horror")
xlabel('time');
ylabel('Value');
title('PSD alpha - P7/P8');