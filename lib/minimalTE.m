function [TEmin, DEL, del] = minimalTE(Gmax, Smax, bmax, Nx, PF, GRAPPA, t_RF90, t_RF180, t_ADCstart, esp)
% Calculate the minimal TE using the widest pulse
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

% pulse width is constrained by the TE, ADC duration, 180 RF pulse, and
% rise time
% delta = @(TE) TE/2 - t_ADC - t_RF180/2 - t_rise;

% inter-pulse duration is constrained by the TE, 90 RF pulse, and 180 RF
% pulse
% Delta = @(TE) TE/2 - t_RF90/2 + t_RF180/2;

% gmax = gamma * Gmax
gmax = Gmax/40*0.0107;

% Simple calculation of b-value without the consideration of rise time
% b = @(TE) gmax^2 * delta(TE).^2 * (Delta(TE)-delta(TE)/3) - bmax;

% Full calculation of b-value with the consideration of rise time. However,
% the triangular pulse is not considered.
% bval = @(TE) gmax^2 * ( 1/30*t_rise^3 - 1/6*t_rise^2*delta(TE) +...
%     t_rise*delta(TE).^2 + 2/3*delta(TE).^3 +...
%     delta(TE).^2.*(Delta(TE)-delta(TE)-t_rise) ) - bmax;

% Calculate the minimal TE
TE0 = (2*t_rise + t_RF180/2 + t_ADCstart) + (0:10:400);        % initialization
TEmin = Inf;
% Full solution of b-value with the consideration of rise time and
% triangular pulse.
fun = @(TE) bval_offset(gmax, TE, t_ADC, t_RF90, t_RF180, t_rise, bmax);
for i = 1:numel(TE0)
    % Calculate minimal TE such that b = bmax
    TEi = fzero(fun, TE0(i));
    DELi = Delta(TEi, t_RF90, t_RF180);
    deli = delta(TEi, t_ADC, t_RF180, t_rise);
    if DELi>0 && deli>0 && TEi<TEmin
        TEmin = TEi;
    end
end
DEL = Delta(TEmin, t_RF90, t_RF180);
del = delta(TEmin, t_ADC, t_RF180, t_rise);

end

function [del, del0] = delta(TE, t_ADC, t_RF180, t_rise)
% pulse width with the consideration of triangular pulse
del0 = TE/2 - t_ADC - t_RF180/2 - t_rise;
if del0 < t_rise
    del = (TE/2 - t_ADC - t_RF180/2)/2;
else
    del = del0;
end

end

function DEL = Delta(TE, t_RF90, t_RF180)
% inter-pulse duration
DEL = TE/2 - t_RF90/2 + t_RF180/2;
end

function b = bval_offset(gmax, TE, t_ADC, t_RF90, t_RF180, t_rise, bmax)
% Full solution of b-value with the consideration of rise time and
% triangular pulse.

[del, del0] = delta(TE, t_ADC, t_RF180, t_rise);
DEL = Delta(TE, t_RF90, t_RF180);

if del0 < t_rise
    gmax = gmax * (del/t_rise);
    t_rise = del;
end

b = gmax^2 * ( 1/30*t_rise^3 - 1/6*t_rise^2*del +...
    t_rise*del.^2 + 2/3*del.^3 +...
    del.^2.*(DEL-del-t_rise) ) - bmax;
end



