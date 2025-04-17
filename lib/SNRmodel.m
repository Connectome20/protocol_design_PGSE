function SNR = SNRmodel(f, T2, TE, Nx, PF, GRAPPA, esp)
% f: volume fraction in multiple compartments
% T2: T2 values in multiple compartments, ms
% TE: echo time, ms
% Nx: # kx in each line
% PF: partial fourier factor
% GRAPPA: GRAPPA acceleration factor
% esp: echo spacing, ms
%
% SNR: signal-to-noise ratio
%
% Author: Hong-Hsi Lee, HLEE84@mgh.harvard.edu

Ny = Nx*PF/GRAPPA;
BW = 1000/esp*Ny;
SNR = f.*exp(-TE./T2)*sqrt(Nx*Ny)/sqrt(BW);
SNR = sum(SNR(:));

end