function DASM = dasm_signal(signal1, signal2, epoch_length, Fs)
    [signal1, ~, ~, ~] = signal_decompose(signal1, Fs, 'alpha');
    [signal2, ~, ~, ~] = signal_decompose(signal2, Fs, 'alpha');
    DASM = differential_entropy(signal1,epoch_length,Fs) - differential_entropy(signal2,epoch_length,Fs);
end