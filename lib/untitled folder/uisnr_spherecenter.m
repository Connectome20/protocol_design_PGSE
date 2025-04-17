function snr_ult = uisnr_spherecenter(fieldstrength,epsilon_rel,sigma,sphereradius)
%UISNR_SPHERECENTER   Ultimate intrinsic SNR at the sphere center.
%   UISNR_SPHERECENTER(fieldstrength,epsilon_rel,sigma,sphereradius) is the
%   ultimate instrinsic SNR at the center of a dielectric sphere in radius
%   (sphereradius), with the relative permittivity (epsilon_rel) and the
%   electric conductivity (sigma) at the magnetic field (filedstrength).
%
%   Input:
%     fieldstrength: field strength (Tesla)
%     epsilon_rel: relative permittivity of imaged body
%         [NaN --> select brain tissue values based on the field strength]
%     sigma: conductivity of imaged body (Siemens/m) 
%         [NaN --> select brain tissue values based on the field strength]
%     sphereradius: radius of dielectric sphere (m)
%   Output:
%     ultimate intrinsic SNR
%
%   (c) Hong-Hsi Lee, Jan 26, 2016

% expsnr_num: experimental SNR scaling, we do not discuss this issue here.
expsnr_num = 1;

% -- determine wavevectors -- 
mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]
c = 3e8;                % speed of light [m]/[s]
epsilon_0 = 1/(mu*c^2); % permittivity [C][V^-1] = [s^4][A^2][m^-2][kg^-2]
omega = 2*pi*42.576e6*fieldstrength; % Larmor frequency [Hz]

% Brain Tissue Properties (Gabriel S, et al., Phys Med Biol, 1996)
fieldset = [1 3 5 7 9 11];
epsilon_rel_brain = [102.5 63.1 55.3 52 50 48.8];
sigma_brain = [0.36 0.46 0.51 0.55 0.59 0.62];
if isnan(epsilon_rel)
    epsilon_rel = spline(fieldset,epsilon_rel_brain,fieldstrength);
else
end

if isnan(sigma)
    sigma = spline(fieldset,sigma_brain,fieldstrength);
else
end
epsilon = epsilon_rel*epsilon_0;

% disp('-------------------------------------------');
% disp(['B_o = ' num2str(fieldstrength) ' [T]']);
% disp(['omega = ' num2str(omega/1E6) ' [MHz]']);
% disp(['sigma = ' num2str(sigma) ' [ohm^-1][m^-1]']);
% disp(['epsilon = ' num2str(epsilon) ' [C][V^-1]']);
% disp('-------------------------------------------');

k_in_squared = omega*mu*(omega*epsilon+1j*sigma);
k_in = sqrt(k_in_squared);

%--------------------------------------------
%        Ultimate intrinsic SNR scaling
%--------------------------------------------

Nproton = 6.691e28;     % number of protons per unit volume in water
gyromag = 2.68e8;       % gyromagnetic ratio for protons [rad/T/sec]
h = 6.626e-34;          % Planck's constant [J sec]
hbar = h/2/pi;
Ispin = 1/2;
k_B = 1.3806503e-23;    % Boltzmann's constant
T = 310;
M_0 = Nproton*((gyromag*hbar)^2)*Ispin*(Ispin+1)*fieldstrength/(3*k_B*T);

usnr_num = omega*M_0/sqrt(4*k_B*T);

% ** experimental SNR scaling is defined in the Appendix, we do not
% discuss this issue here.
snr_num = usnr_num*expsnr_num;

snr_denom = sqrt( 3*pi*sigma*omega^2/abs(k_in)^2*computePsiM_center(k_in,sphereradius) );

snr_ult = snr_num./snr_denom;

end