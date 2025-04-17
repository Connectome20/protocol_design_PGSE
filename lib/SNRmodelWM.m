function SNR = SNRmodelWM(f, T2, TE, Nx, PF, GRAPPA, esp, b, Da)
% f: intra-cellular volume fraction
% T2: T2 values, ms
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
DWI = sqrt(pi/4/b/Da)*erf(sqrt(b*Da));
SNR = f*exp(-TE/T2) * DWI *sqrt(Nx*Ny)/sqrt(BW);
SNR = sum(SNR(:));

end