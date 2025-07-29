function [amp_fts, phs_fts, omega] = fft_forward(dt, signal)
    Fs = 1/dt;        % Sampling frequency (Hz)
    L = length(signal);
    fts = fft(signal) / L;  % Normalize FFT by length
    
    if mod(L,2) == 0
       N_half = L/2 + 1;
    else
       N_half = (L+1)/2;
    end
    
    amp_fts = abs(fts(1:N_half));
    amp_fts(2:end) = amp_fts(2:end) * 2;

    phs_fts = angle(fts(1:N_half))+pi/2;

    omega = (0:N_half-1) * (Fs/L);
end