clc; clear all; close all;

Fs = 128;
addpath 'C:\Users\david\OneDrive\Documentos\GitHub\MindSpark\functions\sigprocfunc';
% Get the mean signal for relax state
subjects = dir('C:\Users\david\OneDrive\Documentos\GitHub\MindSpark\GAMEEMO');

epoch_length = 30; 
features = zeros(1,15);

for i=1:length(subjects)
    if contains(subjects(i).name,')')
        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G2AllChannels.mat');
        str = {'C:\Users\david\OneDrive\Documentos\GitHub\MindSpark\GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);
        calm_AF3 = AF3';
        calm_AF4 = AF4';
        calm_P7 = P7';
        calm_P8 = P8';

        file_str = strcat(replace(subjects(i).name, {'(', ')'}, ''), 'G3AllChannels.mat');
        str = {'C:\Users\david\OneDrive\Documentos\GitHub\MindSpark\GAMEEMO', subjects(i).name, 'Preprocessed EEG Data', '.mat format', file_str};
        str = strjoin(str, '\');
        load(str);

        horror_AF3 = AF3';
        horror_AF4 = AF4';
        horror_P7 = P7';
        horror_P8 = P8';
        
        % low_cf = 1;
        % high_cf = 70;
        % filter_order = 4;
        % [b, a] = butter(filter_order, [low_cf high_cf]/(Fs/2));
        % 
        % calm_AF3 = filtfilt(b,a,calm_AF3);
        % horror_AF3 = filtfilt(b,a,horror_AF3);
        
        eeg_signal_AF3 = [calm_AF3;horror_AF3];
        eeg_signal_AF4 = [calm_AF4;horror_AF4];
        eeg_signal_P7 = [calm_P7; horror_P7];
        eeg_signal_P8 = [calm_P8; horror_P8];

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

        alpha_signalAF = (alpha_signal_AF4 + alpha_signal_AF3) ./2;
        beta_signalAF = (beta_signal_AF4 + beta_signal_AF3) ./2;
        gamma_signalAF = (gamma_signal_AF4 + gamma_signal_AF3) ./2;
        theta_signalAF = (theta_signal_AF4 + theta_signal_AF3) ./2;
        
        
        alpha_signal_P7 = eegfilt(eeg_signal_P7, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_P7 = eegfilt(eeg_signal_P7, Fs, beta_freq(1), beta_freq(2));
        gamma_signal_P7 = eegfilt(eeg_signal_P7, Fs, gamma_freq(1), gamma_freq(2));
        theta_signal_P7 = eegfilt(eeg_signal_P7, Fs, theta_freq(1), theta_freq(2));

        alpha_signal_P8= eegfilt(eeg_signal_P8, Fs, alpha_freq(1), alpha_freq(2));
        beta_signal_P8 = eegfilt(eeg_signal_P8, Fs, beta_freq(1), beta_freq(2));
        gamma_signal_P8 = eegfilt(eeg_signal_P8, Fs, gamma_freq(1), gamma_freq(2));
        theta_signal_P8 = eegfilt(eeg_signal_P8, Fs, theta_freq(1), theta_freq(2));

        alpha_signalP = (alpha_signal_P7 + alpha_signal_P8) ./2;
        beta_signalP = (beta_signal_P7 + beta_signal_P8) ./2;
        gamma_signalP = (gamma_signal_P7 + gamma_signal_P8) ./2;
        theta_signalP = (theta_signal_P7 + theta_signal_P8) ./2;

        for m=1:size(eeg_signal_AF3)
            % AF features
            % Energy alfa/Gamma
            EalphaGamma = signal_energy(alpha_signalAF(m,:), epoch_length, Fs)./signal_energy(gamma_signalAF(m,:), epoch_length, Fs);
            % Energy alfa/beta_ASM
            EalphaBetaASM = signal_energy(alpha_signalAF(m,:), epoch_length, Fs)./asm_signal(beta_signal_AF3,beta_signal_AF4,epoch_length, Fs);
            % Entropy Energy beta
            EEbeta = energy_entropy(beta_signalAF(m,:), epoch_length, Fs);
            % Entropy Energy gamma
            EEgamma = energy_entropy(gamma_signalAF(m,:), epoch_length, Fs);
            % Entropy Energy theta
            EEtheta = energy_entropy(theta_signalAF(m,:), epoch_length,Fs);
            % DASM alpha
            DASMalpha = dasm_signal(alpha_signal_AF3(m,:), alpha_signal_AF4(m,:), epoch_length, Fs);
            % PSD beta
            PSD_beta_AF= psd_signal(beta_signalAF(m,:), epoch_length, Fs, beta_freq); 

            % TP (based on P7 and P8 instead of TP9 and TP10)
            % Energy alfa/gamma
            EalphaGammaP = signal_energy(alpha_signalP(m,:), epoch_length, Fs)./signal_energy(gamma_signalP(m,:), epoch_length, Fs);
            % Energy alfa/beta_ASM
            EalphaBetaASMP = signal_energy(alpha_signalP(m,:), epoch_length, Fs)./asm_signal(beta_signal_P7,beta_signal_P8,epoch_length, Fs);
            % Entropy Energy beta
            EEbetaP = energy_entropy(beta_signalP(m,:), epoch_length, Fs);
            % Entropy Energy gamma
            EEgammaP = energy_entropy(gamma_signalP(m,:), epoch_length,Fs);
            % Entropy Energy theta
            EEthetaP = energy_entropy(theta_signalP(m,:), epoch_length,Fs);
            % DASM alpha
            DASMalphaP = dasm_signal(alpha_signal_P7(m,:), alpha_signal_P8(m,:), epoch_length, Fs);
            % PSD beta
            PSD_beta_P= psd_signal(beta_signalAF(m,:), epoch_length, Fs, beta_freq); 


            samples = [EalphaGamma', EalphaBetaASM', EEbeta', EEgamma', EEtheta', DASMalpha', PSD_beta_AF'];
            samples(:,end+1:end+8) = [EalphaGammaP', EalphaBetaASMP', EEbetaP', EEgammaP', EEthetaP', DASMalphaP', PSD_beta_P', ones(length(PSD_beta_P),1)*(m-1)];
            features(end+1:end+length(PSD_beta_P),:) = samples;
        end
    end
end

features = features(2:end,:);
features = array2table(features);
features.Properties.VariableNames(1:7) = {'EalphaGamma','EalphaBetaASM','EEbeta','EEgamma','EEtheta','DASMalpha','psdBeta'};
features.Properties.VariableNames(8:15) = {'EalphaGammaP','EalphaBetaASMP','EEbetaP','EEgammaP','EEthetaP','DASMalphaP','psdBetaP','Class'};

writetable(features,'features.csv');

%%
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