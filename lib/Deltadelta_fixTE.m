function [DEL, del, G] = Deltadelta_fixTE(Gmax, Smax, Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp, TE, N, M, del_min)
% Calculate the (Delta, delta, Gmax) at a given echo time
% Gmax: maximal gradient strength, mT/m
% Smax: maximal slew rate, T/m/s
% Nx: # kx in each line
% PF: partial fourier factor
% GRAPPA: GRAPPA acceleration factor
% t_RF90: time width of the 90 degree RF pulse, ms
% t_RF180: time width of the 180 degree RF pulse, ms
% t_ADCstart: time of traveling from center k-space to the start of EPI, ms
% esp: echo spacing of EPI, ms
% TE: echo time, ms
% N: # delta for each Delta
% M: # Delta for the widest delta
%
% DEL: Delta, inter-pulse duration, ms
% del: delta, pulse width, ms
% G: maximum of allowable gradient strength, mT/m
%
% Author: Hong-Hsi Lee, HLEE84@mgh.harvard.edu

% Time from the start of EPI to the echo of EPI
t_ADC = Nx*(PF-1/2)/GRAPPA * esp + t_ADCstart;

% Rise time of the diffusion pulse
t_rise = Gmax/Smax;

% time between RF 180 and ADC duration
t2 = TE/2 - t_ADC - t_RF180/2;
if t2 < 0
    error('Echo time <%.2f ms is too short to accomodate the 180 RF pulse and ADC duration.\n', 2*t_ADC + t_RF180);
end

% calculate delta, t_rise, and G
if t2 < 2*t_rise
    % triangular pulse
    if del_min > t2/2
        % set del_min to 0 if it is too large
        del_min = 0;
    end
    Gmax_min = Gmax*(del_min/t_rise);
    dels = linspace(del_min,t2/2,N);
    dels = dels(:);
    t_rises = dels;
    Gmaxs = linspace(Gmax_min,Gmax*(t2/(2*t_rise)),N);
    Gmaxs = Gmaxs(:);
else
    % trapezoid pulse
    dels = linspace(del_min,t2-t_rise,N);
    dels = dels(:);
    t_rises = repmat(t_rise,N,1);
    Gmaxs = repmat(Gmax,N,1);
end

% Delta bounds at the shortest delta
DEL_ub_1 = TE - t_ADC - dels(1) - t_rises(1) - t_RF90/2;
DEL_lb_1 = dels(1) + t_rises(1) + t_RF180;

% Delta bounds at the widest delta
DEL_ub_N = TE - t_ADC - dels(N) - t_rises(N) - t_RF90/2;
DEL_lb_N = dels(N) + t_rises(N) + t_RF180;

% Delta increment in the middle time range
dt2 = (DEL_ub_N - DEL_lb_N)/(M-1);

% Delta increment in the upper and lower time range
dt1 = (DEL_ub_1 - DEL_ub_N)/ceil((DEL_ub_1 - DEL_ub_N)/dt2);

% calculate Delta and delta
DEL = [];
del = [];
G   = [];
for i = 1:N
    deli = dels(i);
    t_risei = t_rises(i);
    Gmaxi = Gmaxs(i);
    DEL_ub = TE - t_ADC - deli - t_risei - t_RF90/2;
    DEL_lb = deli + t_risei + t_RF180;
    DEL1 = [DEL_ub_N     :  dt1 : DEL_ub, DEL_ub]; DEL1 = unique(DEL1);
    DEL2 =  DEL_lb_N+dt2 :  dt2 : DEL_ub_N-dt2;
    DEL3 = [DEL_lb_N     : -dt1 : DEL_lb, DEL_lb]; DEL3 = unique(DEL3);
    DELi = [DEL1, DEL2, DEL3].';
    DELi = sort(DELi);
    DEL  = cat(1, DEL, DELi);
    del  = cat(1, del, deli*ones(numel(DELi),1));
    G    = cat(1, G, Gmaxi*ones(numel(DELi),1));

    % for j = 1:N
    %     k = k+1;
    %     DELi = DELs(j);
    %     del(k) = deli;
    %     DEL(k) = DELi;
    %     G(k)   = Gmaxi;
    % end
end
    
end




