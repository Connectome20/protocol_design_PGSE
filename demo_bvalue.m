%% Demonstration of b-value calculation with the consideration of rise time
close all
clc

G = 500;            % gradient strength, mT/m
S = 600;            % slew rate, T/m/s
g = G/40*0.0107;    % Larmor gradient strength, 1/um/ms
t_rise = G/S;       % rise time, ms
delta = 10;         % pulse width, ms
Delta = 20;         % inter-pulse duration, ms

% triangular pulse
if delta < t_rise
    g = g * (delta/t_rise);
    t_rise = delta;
end

% analytical calculation of b-value
bval_analytical = g^2 * ( 1/30*t_rise^3 - 1/6*t_rise^2*delta +...
    t_rise*delta.^2 + 2/3*delta.^3 +...
    delta.^2.*(Delta-delta-t_rise) );

% create the pulse gradient sequence
Nt = 1000;                      % # time points
dt = (Delta+delta+t_rise)/Nt;   % time step
tt = (1:Nt)*dt;                 % time
gt = zeros(Nt,1);               % gradient sequence
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

% q = integration of g * dt
qt = cumsum(gt)*dt;

% numerical calculation of b-value = integration of q^2 * dt
bval_numerical = sum(qt.^2)*dt;

figure('unit','inch','position',[0 0 10 5]);
subplot(121);
plot(tt, gt);
xlim([0 Delta + delta + t_rise]);
xlabel('time, ms','interpreter','latex','fontsize',20);
ylabel('Larmor gradient strength, 1/um/ms','interpreter','latex','fontsize',20);

subplot(122);
plot(tt, qt);
xlim([0 Delta + delta + t_rise]);
xlabel('time, ms','interpreter','latex','fontsize',20);
ylabel('q-vector, 1/um','interpreter','latex','fontsize',20);
title(sprintf('$b_{\\rm analytical}$=%.4f, $b_{\\rm numerical}$=%.4f',...
    bval_analytical, bval_numerical),'interpreter','latex','fontsize',20);




