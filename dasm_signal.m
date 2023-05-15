function DASM = dasm_signal(signal1, signal2, epoch_length, Fs)
    DASM = differential_entropy(signal1,epoch_length,Fs) - differential_entropy(signal2,epoch_length,Fs);
end