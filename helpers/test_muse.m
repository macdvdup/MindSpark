% *Test Muse Headband* 
%
%     Developers: David Machado, Érica Gomes, José Almeida, Raquel Sousa
%     Neural Engineering Course, FEUP 2023
%
%     This code allows to connect and save data from the Muse Headband. It
%     was used to test different EEG preprocessing methods, features and
%     classification models.

clear all; close all;

% Define the acquisition time
acquisition_time = 65;

% Connect to Muse Headband
mu = MuseUdp();

% Initialize the variables
cont = true;
vector_tp9 = [];
vector_af7 = [];
vector_af8 = [];
vector_tp10 = [];
vector_raux = [];

% Acquistion loop using tic and toc methods from matlab
tic
while cont
    % Tries to get a eeg sample
    [data, timestamp, success] = mu.get_eeg_sample();

    % If it was able to get data, then store it
    if success
        vector_tp9(end+1) = data(1);
        vector_af7(end+1) = data(2);
        vector_af8(end+1) = data(3);
        vector_tp10(end+1) = data(4);
        vector_raux(end+1) = data(5);
    end
    
    if toc > acquisition_time
        cont = false;
    end
end

% Stores the data into a .csv file
M = [vector_tp9, vector_af7, vector_af8, vector_tp10, vector_raux];
writematrix(M, 'data_andreia_2.csv')