%% Calculation of minimal TE and the SNR benefit for a given delta
close all
filePath = matlab.desktop.editor.getActiveFilename;
filePath = fileparts(filePath);
addpath(genpath(fullfile(filePath,'..','lib')))

% load experimental data
file = load(fullfile(filePath,'snr.mat'));

Gmax = [480 300 80];        % maximal gradient strength, mT/m
Gmax_lg = [500 300 80];     % maximal gradient strength in legend
Smax = [600 80 80];         % maximal slew rate, T/m/s
bmax = 0.1:0.1:40;          % maximal b-value, ms/um2
delta = 8;                  % pulse width, ms
Nx = 110;                   % # ky lines
PF = 6/8;                   % partial fourier factor
GRAPPA = 2;                 % GRAPPA acceleration factor
t_RF90 = 2.1;               % time width of the 90 degree RF pulse, ms
t_RF180 = 3.4;              % time width of the 180 degree RF pulse, ms
t_ADCstart = 0.2;           % time of traveling from center k-space to the start of EPI, ms
t_DeadTime = 0.3;           % dead time after 2nd diffusion gradient pulse
t_ADCstart = t_ADCstart +...% include the dead time into imaging prewinder
             t_DeadTime;
esp = 0.4;                  % echo spacing of EPI, ms

% Time of reference scan for N/2 ghost correction right after 90 degree RF pulse, ms
% Typical values: prewinder = 0.34 ms, 3 reference lines = 3*esp, rewinder = 0.22 ms
t_REF = 0.34+3*esp+0.22;
% Include the time width of reference scan into time width of 90 degree RF pulse
t_RF90 = t_RF90 + 2*t_REF;

% calculate minimal TE at fixed delta
Nbmax = numel(bmax);
NGmax = numel(Gmax);
TEmin = zeros(NGmax,Nbmax);
DEL   = zeros(NGmax,Nbmax);
del   = zeros(NGmax,Nbmax);
for i = 1:NGmax
    for j = 1:Nbmax
        [TEmin(i,j), DEL(i,j), del(i,j)] = minimalTE_fixdelta(Gmax(i), Smax(i), bmax(j), Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp, delta);
    end
end

% plot b-value vs minimal TE
figure('unit','inch','position',[0 0 10 5]);
subplot(121);
cmap = colormap('lines');
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,TEmin(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax_lg(i),round(Smax(i)));
end
plot(5, 176,'.','Color',cmap(3,:),'markersize',10);
plot([5 10 20 30 40],[46 55 79 106 128]-6,'.','Color',cmap(2,:),'markersize',10)
plot([5 10 20 30 40],[40 40 48 57 65]-6,'.','Color',cmap(1,:),'markersize',10)

xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('TE$_{\rm min}$','interpreter','latex','fontsize',20);
ylim([0 200]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','north','box','off');


% calculate theoretical SNR based on T2 value in white matter
f  = 1;                         % volume fraction
B0_exp = [3 4 7];               % B0 field strength, T
T2_exp = [80 58 46];            % T2 values, ms, Wansapura et al. JMRI 1999; Yacoub et al. MRM 2003
T2_of_B0 = @(x) 1./interp1(1./B0_exp, 1./T2_exp, 1./x, 'spline');
SNR = zeros(NGmax,Nbmax);
for i = 1:NGmax
    for j = 1:Nbmax
        B0 = 3;
        T2 = T2_of_B0(B0);
        SNR(i,j) = SNRmodel(f, T2, TEmin(i,j), Nx, PF, GRAPPA, esp);
    end
end

% plot b(TE_min) vs SNR gain wrt C1
subplot(122);
hold on;
clear h lgtxt
cmap = colormap('lines');
for i = 1:NGmax
    h(i) = plot(bmax, SNR(i,:)./SNR(2,:),'-','linewidth',1,'color',cmap(i,:));
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax_lg(i),round(Smax(i)));
end

% C1 result
[b_C1, I1] = sort(file.b_C1); snr_b0_C1 = file.snr_b0_C1(I1);

% C2 result
[b_C2, I2] = sort(file.b_C2); snr_b0_C2 = file.snr_b0_C2(I2);

% plot SNR ratio
plot(b_C2, snr_b0_C2./snr_b0_C1, '.', 'markersize', 10, 'color', cmap(1,:));
plot(b_C1, snr_b0_C1./snr_b0_C1, '.', 'markersize', 10, 'color', cmap(2,:));
plot(file.b_prisma, file.snr_b0_prisma/snr_b0_C1(1), '.', 'markersize', 10, 'color', cmap(3,:));

xlabel('$b$(TE$_{\rm min}$)','interpreter','latex','fontsize',20);
ylabel('SNR grain wrt C1','interpreter','latex','fontsize',20);
ylim([0 2.5]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt,...
    'interpreter','latex','fontsize',12,'location','north','box','off')

%% SNR comparison of 3T vs 7T scanner
% calculate theoretical SNR based on T2 value in white matter
f  = 1;                         % volume fraction
B0_exp = [3 4 7];               % B0 field strength, T
T2_exp = [80 58 46];            % T2 values, ms, Wansapura et al. JMRI 1999; Yacoub et al. MRM 2003
T2_of_B0 = @(x) 1./interp1(1./B0_exp, 1./T2_exp, 1./x, 'spline');

B0 = [3 7];
NB0 = numel(B0);
SNR = zeros(NGmax,Nbmax,NB0);
for i = 1:NGmax
    for j = 1:Nbmax
        for k = 1:NB0
            % calculate ultimate intrinsic SNR at each B0 field strength
            sphereradius = 0.07; % head radius, meter
            snr_ult = uisnr_spherecenter(B0(k), NaN, NaN, sphereradius);

            T2 = T2_of_B0(B0(k));
            SNR(i,j,k) = snr_ult*SNRmodel(f, T2, TEmin(i,j), Nx, PF, GRAPPA, esp);
        end
    end
end

figure('unit','inch','position',[0 0 15 5]);
% plot b-value vs minimal TE
subplot(131);
cmap = colormap('lines');
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,TEmin(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax_lg(i),round(Smax(i)));
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
        lgtxt{l} = sprintf('%u mT/m, %u T/m/s, %u T',Gmax_lg(i),round(Smax(i)),B0(k));
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



