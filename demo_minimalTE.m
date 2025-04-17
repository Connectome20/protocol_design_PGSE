addpath(genpath('/Users/honghsilee/Documents/GitHub/C2_protocoldesign/lib'));

%% Calculation of minimal TE and the SNR benefit
% close all

Gmax = [80 300 500];        % maximal gradient strength, mT/m
Smax = [80 80 600];         % maximal slew rate, T/m/s
bmax = 0.1:0.1:40;          % maximal b-value, ms/um2
Nx = 110;                   % # kx in each line
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

Nbmax = numel(bmax);
NGmax = numel(Gmax);
TEmin = zeros(NGmax,Nbmax);
DEL   = zeros(NGmax,Nbmax);
del   = zeros(NGmax,Nbmax);
for i = 1:NGmax
    for j = 1:Nbmax
        [TEmin(i,j), DEL(i,j), del(i,j)] = minimalTE(Gmax(i), Smax(i), bmax(j), Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp);
    end
end

figure('unit','inch','position',[0 0 10 10]);
subplot(221);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,TEmin(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('TE$_{\rm min}$','interpreter','latex','fontsize',20);
ylim([20 120]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','northwest','box','off');

subplot(223);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,DEL(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('$\Delta$','interpreter','latex','fontsize',20);
ylim([0 60]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','northwest','box','off');

subplot(224);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,del(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('$\delta$','interpreter','latex','fontsize',20);
ylim([0 60]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','northwest','box','off');

f = 0.7;                % intra-cellular volume fraction
T2 = 90;                % intra-cellular T2, ms
SNR = zeros(NGmax,Nbmax);
for i = 1:NGmax
    for j = 1:Nbmax
        SNR(i,j) = SNRmodel(f, T2, TEmin(i,j), Nx, PF, GRAPPA, esp);
    end
end

subplot(222);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,SNR(i,:)./SNR(1,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('relative SNR wrt 40 mT/m','interpreter','latex','fontsize',20);
ylim([1 2.2]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','northwest','box','off');
















%%
% za = 1.64;
% n = 64;
% SNR = 50;
% sigmabar = za/SNR/sqrt(n);
% 
% D0 = 2;
% delta = 10;
% G = 500;
% g = G/40*0.0107;
% Delta = 20;
% Da = 1.7;
% 
% b = g^2 * delta^2 * (Delta-delta/3);
% hA = sqrt(pi/4/b/Da)*erf(sqrt(b*Da));
% 
% d_min = (768/7 * sigmabar * D0 /delta /g^2)^(1/4) * hA^(-1/4);

