function ASM = asm_signal(signal1, signal2, epoch_length, Fs)
    ASM = signal_energy(signal1,epoch_length,Fs) - signal_energy(signal2,epoch_length,Fs);
end