function [t] = preProcessEEG(vector, cutOffFreq, Fs)
    filter_order = 6;
    [b, a] = butter(filter_order, [cutOffFreq(1) cutOffFreq(2)]/(Fs/2));
    t = filtfilt(b,a,vector);
end