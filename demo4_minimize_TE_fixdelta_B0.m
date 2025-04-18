%% Demo 4: Calculation of the shortest TE of PGSE for a fixed pulse width at 3T vs 7T
% Reference: Ramos-Llorden, Lee, ..., Huang, Nature BME, 2025
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);
addpath(genpath(fullfile(root,'lib')));

% Diffusion spec for Connectome 2.0, Connectome 1.0, and clinical scanners
Gmax = [500 300 80];        % maximal gradient strength, mT/m
Smax = [600 80  80];        % maximal slew rate, T/m/s
bmax = 0.1:0.1:40;          % maximal b-value, ms/um2
delta = 8;                  % pulse width, ms

% Spin echo sequence parameters
Nx = 110;                   % image matrix size, Nx by Nx
PF = 6/8;                   % partial fourier factor
GRAPPA = 2;                 % GRAPPA acceleration factor
t_RF90 = 2.1;               % time width of the 90 degree excitation RF pulse, ms
t_RF180 = 3.4;              % time width of the 180 degree refocusing RF pulse, ms
t_ADCstart = 0.2;           % time to travel from center k-space to the start of EPI (imaging prewinder), ms
t_DeadTime = 0.3;           % dead time after the 2nd diffusion gradient pulse, ms
t_ADCstart = t_ADCstart +...% Here we include the dead time into imaging prewinder.
             t_DeadTime;
esp = 0.4;                  % echo spacing of EPI, ms

% Time of reference scan for N/2 ghost correction right after 90 degree excitation RF pulse, ms
% Typical values for reference scan: prewinder = 0.34 ms, 3 reference lines = 3*esp, rewinder = 0.22 ms
t_REF = 0.34 + 3*esp + 0.22;

% Here we include the time width of reference scan into time width of 90 degree excitation RF pulse.
t_RF90 = t_RF90 + 2*t_REF;

% Calculate the shortest TE for a fixed pulse width
Nbmax = numel(bmax);        % # b-value
NGmax = numel(Gmax);        % # gradient strength
TEmin = zeros(NGmax,Nbmax); % shortest echo time, ms
DEL   = zeros(NGmax,Nbmax); % inter-pulse duration (diffusion time), ms
del   = zeros(NGmax,Nbmax); % pulse width, ms
for i = 1:NGmax
    for j = 1:Nbmax
        [TEmin(i,j), DEL(i,j), del(i,j)] = minimalTE_fixdelta(...
            Gmax(i), Smax(i), bmax(j),...
            Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp, delta);
    end
end

% Calculate SNR based on T2 value in white matter at 3T vs 7T
f  = 1;                         % volume fraction of multiple compartments, sum(f) = 1
B0_exp = [3 4 7];               % B0 field strength, T
T2_exp = [80 58 46];            % T2 values, ms (Wansapura et al. JMRI 1999; Yacoub et al. MRM 2003)
T2_of_B0 = @(x) 1./interp1(1./B0_exp, 1./T2_exp, 1./x, 'spline');

B0 = [3 7];
NB0 = numel(B0);
SNR = zeros(NGmax,Nbmax,NB0);
for i = 1:NGmax
    for j = 1:Nbmax
        for k = 1:NB0
            % Calculate ultimate intrinsic SNR at a dielectric sphere
            % center at B0 (Lee, et al., MRM 2018)
            sphereradius = 0.07; % head radius, meter
            snr_ult = uisnr_spherecenter(B0(k), NaN, NaN, sphereradius);

            T2 = T2_of_B0(B0(k));
            SNR(i,j,k) = snr_ult*SNRmodel(f, T2, TEmin(i,j), Nx, PF, GRAPPA, esp);
        end
    end
end

figure('unit','inch','position',[0 0 15 5]);
% plot b-value vs the shortest TE
subplot(131);
cmap = colormap('lines');
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,TEmin(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),round(Smax(i)));
end
xlabel('$b$, ms/$\mu$m$^2$','interpreter','latex','fontsize',20);
ylabel('TE$_{\rm min}$, ms','interpreter','latex','fontsize',20);
ylim([0 200]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','north','box','off');

% plot B0 vs T2 in WM
subplot(132);
cmap = colormap('flag');
hold on;
clear h lgtxt
B0s = 3:0.1:7;
T2s = T2_of_B0(B0s);
plot(B0s, T2s, 'k-', 'linewidth', 1);
for i = 1:numel(B0_exp)
    h(i) = plot(B0_exp(i), T2_exp(i), 'ok', 'MarkerFaceColor', cmap(i,:), 'MarkerSize', 8);
end
lgtxt = {'Wansapura et al., JMRI 1999', 'Yacoub et al., MRM 2003', 'Yacoub et al., MRM 2003'};
legend(h, lgtxt,...
    'interpreter','latex','fontsize',12,'location','north','box','off')
box on;
grid on;
xlim([3 7]);
ylim([40 90]);
yticks(40:10:90)
pbaspect([1 1 1]);
xlabel('$B_0$, T','interpreter','latex','fontsize',20);
ylabel('$T_2$ in white matter, ms','interpreter','latex','fontsize',20);

% plot SNR gain wrt C1, 3T
subplot(133);
hold on;
clear h lgtxt
cmap = colormap('lines');
ls = {'-','--'};
l = 0;
for i = 1:NGmax
    for k = 1:NB0
        l = l+1;
        h(l) = plot(bmax, squeeze(SNR(i,:,k))./SNR(2,:,1),ls{k},'linewidth',1,'color',cmap(i,:));
        lgtxt{l} = sprintf('%u mT/m, %u T/m/s, %u T',Gmax(i),round(Smax(i)),B0(k));
    end
end
xlabel('$b$(TE$_{\rm min}$, ms/$\mu$m$^2$)','interpreter','latex','fontsize',20);
ylabel('SNR gain wrt C1, 3T','interpreter','latex','fontsize',20);
ylim([0 5]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt,...
    'interpreter','latex','fontsize',12,'location','north','box','off')


