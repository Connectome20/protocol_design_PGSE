%% Demo 1: b-value calculation for trapzoidal pulsed-gradient
% Reference: Ramos-Llorden, Lee, ..., Huang, Nature BME, 2025
% Author: Hong-Hsi Lee (0000-0002-3663-6559)

clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);
addpath(genpath(fullfile(root,'lib')));

G = 500;            % gradient strength (G), mT/m
S = 600;            % slew rate (S), T/m/s
g = G/40*0.0107;    % Larmor gradient strength (gamma*G), 1/um/ms
t_rise = G/S;       % rise time, ms
delta = 10;         % pulse width (FWHM of trapzoidal pulse), ms
Delta = 20;         % inter-pulse duration, ms

% When the pulse width is shorter than rise time, the pulse is triangular
if delta < t_rise
    g = g * (delta/t_rise);
    t_rise = delta;
end

% Analytical solution of b-value for trapzoidal pulsed-gradient sequence
bval_analytical = g^2 * ( 1/30*t_rise^3 - 1/6*t_rise^2*delta +...
    t_rise*delta.^2 + 2/3*delta.^3 +...
    delta.^2.*(Delta-delta-t_rise) );

% Create the pulse gradient sequence
Nt = 1000;                      % # time points
dt = (Delta+delta+t_rise)/Nt;   % time step, ms
tt = (1:Nt)*dt;                 % time, ms
gt = zeros(Nt,1);               % Larmor gradient waveform, 1/um/ms
for i = 1:Nt
    ti = i*dt;
    if ti <= t_rise
        gt(i) = ti/t_rise*g;
    elseif ti <= delta
        gt(i) = g;
    elseif ti <= delta+t_rise
        gt(i) = g-(ti-delta)/t_rise*g;
    elseif ti <= Delta
        gt(i) = 0;
    elseif ti <= Delta+t_rise
        gt(i) = -(ti-Delta)/t_rise*g;
    elseif ti <= Delta + delta
        gt(i) = -g;
    elseif ti <= Delta + delta + t_rise
        gt(i) = -( g-(ti-Delta-delta)/t_rise*g );
    end
end

% Diffusion wavevector q(t) = integration of g(t) * dt, 1/um
qt = cumsum(gt)*dt;

% Numerical calculation of b-value = integration of q(t)^2 * dt, ms/um2
bval_numerical = sum(qt.^2)*dt;

figure('unit','inch','position',[0 0 10 5]);
subplot(121);
plot(tt, gt, 'linewidth', 1);
xlim([0 Delta + delta + t_rise]);
xlabel('time, ms','interpreter','latex','fontsize',20);
ylabel('Larmor gradient strength, ($\mu$m$\cdot$ms)$^{-1}$','interpreter','latex','fontsize',20);

subplot(122);
plot(tt, qt, 'linewidth', 1);
xlim([0 Delta + delta + t_rise]);
xlabel('time, ms','interpreter','latex','fontsize',20);
ylabel('$q$-vector, $\mu$m$^{-1}$','interpreter','latex','fontsize',20);
title(sprintf('$b_{\\rm analytical}$=%.4f, $b_{\\rm numerical}$=%.4f',...
    bval_analytical, bval_numerical),'interpreter','latex','fontsize',20);




