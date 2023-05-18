function ASM = asm_signal(signal1, signal2, epoch_length, Fs)
    [signal1, ~, ~, ~] = signal_decompose(signal1, Fs, 'beta');
    [signal2, ~, ~, ~] = signal_decompose(signal2, Fs, 'beta');
    ASM = signal_energy(signal1,epoch_length,Fs) - signal_energy(signal2,epoch_length,Fs);
end