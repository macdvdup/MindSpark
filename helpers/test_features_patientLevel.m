% Test Features at Subject Level
%
%     Developers: David Machado, Érica Gomes, José Almeida, Raquel Sousa
%     Neural Engineering Course, FEUP 2023
%
%     This code allows to see the different in the features obtained between 
%     Horror and Calm Game. With this, we can see if there any major
%     differences that helps us to directly classify the patient state.
%     The output obtained is a set of plots with the patient number (1 to
%     28) in x-axis and the mean value for the respective feature in the
%     y-axis.

clc; clear all; close all;

% Set the Sampling Frequency of the dataset
Fs = 128;

% Add EGGlab functions to path
addpath 'D:\Neuro\Projeto\MindSpark\functions\sigprocfunc'

% Get the mean signal for relax state
subjects = dir('GAMEEMO');

% Set the epoch duration
epoch_length = 30; % in seconds

% Initialize the variables
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

% Determine the mean value of each feature for each patient
for i=1:length(subjects)
    if contains(subjects(i).name,')')

        % Load the calm signal of the subject
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        % Use the AF3 and AF4 information
        calm_AF3 = AF3';
        calm_AF4 = AF4';
        
        % Load the horror signal of the subject
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        horror_AF3 = AF3';
        horror_AF4 = AF4';
        
        % Applies a band pass filter (Butterworth filter of 4th order)
        low_cf = 1;
        high_cf = 50;
        filter_order = 4;
        [b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
        calm_AF3 = filtfilt(b,a,calm_AF3);
        horror_AF3 = filtfilt(b,a,horror_AF3);
        calm_AF4 = filtfilt(b,a,calm_AF4);
        horror_AF4 = filtfilt(b,a,horror_AF4);
        
        % Create a AF3 and AF4 matrix
        eeg_signal_AF3 = [calm_AF3; horror_AF3];
        eeg_signal_AF4 = [calm_AF4; horror_AF4];
        
        % Decompose the signal into alpha, theta, beta and gamma bands
        % using the eegfilt method from EEGLab
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
        
        % Compute the polar channel by avering the AF3 and AF4 signals
        alpha_signal = (alpha_signal_AF4 + alpha_signal_AF3) ./2;
        beta_signal = (beta_signal_AF4 + beta_signal_AF3) ./2;
        gamma_signal = (gamma_signal_AF4 + gamma_signal_AF3) ./2;
        theta_signal = (theta_signal_AF4 + theta_signal_AF3) ./2;
        
        % Compute the features for each channel and state and stores it
        for m=1:size(eeg_signal_AF3)
            EalphaGamma(m,i) = mean(signal_energy(alpha_signal(m,:), epoch_length, Fs)./signal_energy(gamma_signal(m,:), epoch_length, Fs));
            EalphaBetaASM(m,i) = mean(signal_energy(alpha_signal(m,:), epoch_length, Fs)./asm_signal(beta_signal_AF3,beta_signal_AF4,epoch_length, Fs));
            EEbeta(m,i) = mean(energy_entropy(beta_signal(m,:), epoch_length, Fs));
            EEgamma(m,i) = mean(energy_entropy(gamma_signal(m,:), epoch_length, Fs));
            EEtheta(m,i) = mean(energy_entropy(theta_signal(m,:), epoch_length,Fs));
            DASMalpha(m,i) = mean(dasm_signal(alpha_signal_AF3(m,:), alpha_signal_AF4(m,:), epoch_length, Fs));
            psd_patients(m,i) = mean(psd_signal(beta_signal(m,:), epoch_length, Fs, alpha_freq));
        end
    end
end

% Plots the obtained results
t = linspace(1,length(EalphaGamma), length(EalphaGamma));
figure
plot(t, EalphaGamma);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('E alpha/E gamma - P7/P8');

figure
plot(t, EalphaBetaASM);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('E alpha/ASM beta - P7/P8');

figure
plot(t, EEbeta);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('EE beta - P7/P8');

figure
plot(t, EEgamma);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('EE gamma - P7/P8');

figure
plot(t, EEtheta);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('EE theta - P7/P8');

figure
plot(t, DASMalpha);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('DASM alpha - P7/P8');

figure
plot(t, psd_patients);
legend("Calm","Horror")
xlabel('Patient Number');
ylabel('Value');
title('PSD alpha - P7/P8');