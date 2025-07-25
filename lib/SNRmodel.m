function SNR = SNRmodel(f, T2, TE, Nx, PF, GRAPPA, esp)
% f:  volume fraction in multiple compartments, sum(f) = 1
% T2: T2 values in multiple compartments, ms
% TE: echo time, ms
% Nx: image matrix size, Nx by Nx
% PF: partial fourier factor
% GRAPPA: GRAPPA acceleration factor
% esp: echo spacing or ADC duration for each kx line, ms
%
% SNR: signal-to-noise ratio
%
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

Ny = Nx*PF/GRAPPA;
BW = 1000/esp*Nx;
SNR = f.*exp(-TE./T2)*sqrt(Nx*Ny)/sqrt(BW);
SNR = sum(SNR(:));

end