%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the spherical Bessel function. When the argument 
%  is smaller than a specified threshold series expansion is used.
%
%  INPUT:
%  nu: order of the spherical Bessel function
%  z: argument of the spherical Bessel function
%
%  OUTPUT:
%  spherbesselj: Spherical Bessel function
%
%  Name: spherbessJ
%  Author: Riccardo Lattanzi
%  Created: 16 March 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spherbesselj = spherbessJ(nu, z)

threshold = 0.1;

i1_ind = find(abs(z) >= threshold);
i2_ind = find(abs(z) < threshold);

if length(i1_ind) > 0 % standard expression
    z1 = z(i1_ind); y1 = zeros(length(i1_ind),1);
    y1 = sqrt(pi/2) .* sqrt(1./z1) .* besselj((nu + 1/2), z1);
    spherbesselj(i1_ind)=y1;
end
if length(i2_ind) > 0 % series expansion
    z2 = z(i2_ind); y2 = zeros(length(i2_ind),1);
    tmp_1 = gamma(3/2 + nu);
    tmp_2 = sqrt(pi);
    y2 = (z2.^(1/2 + nu)).*( (2.^(-1-nu)).*tmp_2./(tmp_1.*sqrt(z2)) - ...
        (2.^(-3-nu)).*tmp_2.*(z2.^(3/2))./(tmp_1.*(3/2 + nu)) + ...
        (2.^(-6-nu)).*tmp_2.*(z2.^(7/2))./(tmp_1.*(3/2 + nu).*(5/2 + nu)) - ...
        (2.^(-8-nu)).*tmp_2.*(z2.^(11/2))./(tmp_1.*(3/2 + nu).*(5/2 + nu).*(7/2 + nu)) );
    spherbesselj(i2_ind)=y2;
end

% spherbesselj = sqrt(pi/2) .* sqrt(1./z) .* besselj((nu + 1/2), z);

% END spherbessJ.m