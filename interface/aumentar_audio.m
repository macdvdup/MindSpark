% Specify the file paths of the two audio files
file = 'audio.mp3';

[audio, fs] = audioread(file);

combinedAudio = [audio; audio; audio; audio; audio];

outputFile = 'audio_maior.mp3';
audiowrite(outputFile, combinedAudio, fs);

disp('Audio files combined successfully.');