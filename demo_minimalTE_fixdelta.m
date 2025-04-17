%% Calculation of minimal TE and the SNR benefit for a given delta
close all
% addpath('/Users/hl624/Documents/GitHub/C2_protocoldesign/lib');
addpath('/Users/honghsilee/Documents/GitHub/C2_protocoldesign/lib');

Gmax = [80 300 480];        % maximal gradient strength, mT/m
Smax = [62 80 600];         % maximal slew rate, T/m/s
bmax = 0.1:0.1:40;          % maximal b-value, ms/um2
delta = 8;                  % pulse width, ms
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
        [TEmin(i,j), DEL(i,j), del(i,j)] = minimalTE_fixdelta(Gmax(i), Smax(i), bmax(j), Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp, delta);
    end
end

figure('unit','inch','position',[0 0 10 10]);
subplot(221);
cmap = colormap('lines');
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,TEmin(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
plot(5, 176,'.','Color',cmap(1,:),'markersize',10);
plot([5 10 20 30 40],[46 55 79 106 128]-6,'.','Color',cmap(2,:),'markersize',10)
plot([5 10 20 30 40],[40 40 48 57 65]-6,'.','Color',cmap(3,:),'markersize',10)

xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('TE$_{\rm min}$','interpreter','latex','fontsize',20);
ylim([0 200]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','north','box','off');

subplot(223);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,DEL(i,:),'-','linewidth',1,'Color',cmap(i,:));
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
plot(5, 150,'.','Color',cmap(1,:),'markersize',10);
plot([5 10 20 30 40],[17.9 27.2 51.5 76 100],'.','Color',cmap(2,:),'markersize',10)
plot([5 10 20 30 40],[15 15 23 32.2 39.8],'.','Color',cmap(3,:),'markersize',10)

xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('$\Delta$','interpreter','latex','fontsize',20);
ylim([0 100]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','north','box','off');


subplot(224);
hold on;
clear h lgtxt
for i = 1:NGmax
    h(i) = plot(bmax,del(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('$\delta$','interpreter','latex','fontsize',20);
ylim([0 100]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','north','box','off');

f = [0.48, 0.5, 0.02];                % intra-/extra-cellular/CSF volume fraction
T2 = [70, 100, 2000];                  % intra-/extra-cellular/CSF T2, ms
SNR = zeros(NGmax,Nbmax);
for i = 1:NGmax
    for j = 1:Nbmax
        SNR(i,j) = SNRmodel(f, T2, TEmin(i,j), Nx, PF, GRAPPA, esp);
    end
end

subplot(222);
hold on;
clear h lgtxt
SNR0 = 30/SNR(3,1);
for i = 1:NGmax
    h(i) = plot(bmax,SNR0*SNR(i,:),'-','linewidth',1);
    lgtxt{i} = sprintf('%u mT/m, %u T/m/s',Gmax(i),Smax(i));
end
xlabel('$b$','interpreter','latex','fontsize',20);
ylabel('SNR','interpreter','latex','fontsize',20);
ylim([0 35]);
box on;
grid on;
pbaspect([1 1 1]);
legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','south','box','off');
















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

