% Specify the file paths of the two audio files
file1 = 'audio.mp3';

% Read the audio data from the files
[audio1, fs1] = audioread(file1);

% Combine the audio files by concatenating their data
combinedAudio = [audio1; audio1; audio1; audio1; audio1];

% Write the combined audio to a new file
outputFile = 'audio_maior.mp3';
audiowrite(outputFile, combinedAudio, fs1);

disp('Audio files combined successfully.');