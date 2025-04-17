function b = bval_trapezoid(Gmax, Smax, DEL, del)
% Full solution of b-value with the consideration of rise time and
% triangular pulse.
% Gmax: gradient strength, mT/m
% Smax: slew rate, T/m/s
% DEL: inter-pulse duration, ms
% del: pulse width, ms
%
% Author: Hong-Hsi Lee, HLEE84@mgh.harvard.edu

gmax = Gmax/40*0.0107;
t_rise = Gmax./Smax;

b = gmax.^2 .* ( 1/30*t_rise.^3 - 1/6*t_rise.^2.*del +...
    t_rise.*del.^2 + 2/3*del.^3 +...
    del.^2.*(DEL-del-t_rise) );
end
