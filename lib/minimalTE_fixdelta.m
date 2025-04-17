function [TEmin, DEL, del] = minimalTE_fixdelta(Gmax, Smax, bmax, Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp, delta)
% Calculate the minimal echo time at a given pulse width
% Gmax: maximal gradient strength, mT/m
% Smax: maximal slew rate, T/m/s
% bmax: maximal b-value, ms/um2
% Nx: # ky line
% PF: partial fourier factor
% GRAPPA: GRAPPA acceleration factor
% t_RF90: time width of the 90 degree RF pulse, ms
% t_RF180: time width of the 180 degree RF pulse, ms
% t_ADCstart: time of traveling from center k-space to the start of EPI, ms
% esp: echo spacing of EPI, ms
%
% TEmin: minimal echo time, ms
% DEL: Delta, inter-pulse duration, ms
% del: delta, pulse width, ms
%
% Author: Hong-Hsi Lee, HLEE84@mgh.harvard.edu

% Time from the start of EPI to the echo of EPI
t_ADC = Nx*(PF-1/2)/GRAPPA * esp + t_ADCstart;

% Rise time of the diffusion pulse
t_rise = Gmax/Smax;

% gmax = gamma * Gmax
gmax = Gmax/40*0.0107;

% adjust rise time if rise time is longer than delta
if t_rise > delta
    t_rise = delta;
    gmax = gmax*(delta/t_rise);
end

% Lower bound of TE, determined by delta, rise time, and ADC time
TE_lb  = @(G) 2*(t_RF180/2 + delta + G/Smax + t_ADC);

% Calculate the gradient strength for the case of the 1st pulse right after
% RF 90 and the 2nd pulse right after RF 180 and right before ADC
fun_G0 = @(G) bval_offset(G/40*0.0107, G/Smax, bmax, Delta(TE_lb(G), t_RF90, t_RF180), delta);
G0     = fzero(fun_G0, Gmax/2);

if G0 <= Gmax
    % Calculate the shortest possible Delta for the case of the 2nd pulse
    % right after RF180 and right before ADC
    TEmin = TE_lb(G0);
    fun_G1 = @(G) bval_offset(G/40*0.0107, G/Smax, bmax, Delta_lb(delta, G/Smax, t_RF180), delta);
    G1 = fzero(fun_G1, Gmax/2);
    DEL_lb = Delta_lb(delta, G1/Smax, t_RF180);
    if G1<=Gmax
        DEL = DEL_lb;
    else
        fun_DEL = @(DEL) bval_offset(gmax, t_rise, bmax, DEL, delta);
        DEL = fzero(fun_DEL, DEL_lb);
    end
else
    % Calculate the minimal TE
    TE0 = TE_lb(Gmax):20:TE_lb(Gmax)+100;        % initialization
    TEmin = Inf;
    % Full solution of b-value with the consideration of rise time and
    % triangular pulse.
    fun_TE = @(TE) bval_offset(gmax, t_rise, bmax, Delta_ub(TE, t_ADC, t_RF90, t_rise, delta), delta);
    for i = 1:numel(TE0)
        % Calculate minimal TE such that b = bmax
        TEi = fzero(fun_TE, TE0(i));
        DELi = Delta_ub(TEi, t_ADC, t_RF90, t_rise, delta);
        if DELi>0 && TEi<TEmin
            TEmin = TEi;
            DEL = DELi;
        end
    end
    % DEL = Delta_ub(TEmin, t_ADC, t_RF90, t_rise, delta);
end

del = delta;

end

function DEL = Delta(TE, t_RF90, t_RF180)
% inter-pulse duration
DEL = TE/2 - t_RF90/2 + t_RF180/2;
end

function DEL = Delta_ub(TE, t_ADC, t_RF90, t_rise, delta)
% inter-pulse duration
DEL = TE - t_ADC - t_RF90/2 - delta - t_rise;
end

function DEL = Delta_lb(delta, t_rise, t_RF180)
% lower bound of inter-pulse duration
DEL = delta + t_rise + t_RF180;
end

function b = bval_offset(gmax, t_rise, bmax, DEL, del)
% Full solution of b-value with the consideration of rise time and
% triangular pulse.

b = gmax^2 * ( 1/30*t_rise^3 - 1/6*t_rise^2*del +...
    t_rise*del.^2 + 2/3*del.^3 +...
    del.^2.*(DEL-del-t_rise) ) - bmax;
end




