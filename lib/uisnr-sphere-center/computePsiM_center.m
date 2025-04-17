%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the value of the elements of Psi M using exact
%  solution.
%
%  INPUT:
%  k_in: k_in^2 = omega^2*epsilon*mu + 1j*omega*mu_0*sigma, mu ~ mu_0
%  a: radius of the spherical phantom
%
%  OUTPUT:
%  PsiMvalue: element of Psi M corresponding to l = 1
%
%  Name: computePsiM_center
%  Author: Hong Hsi Lee
%  Created: Jan 25, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PsiMvalue = computePsiM_center(k_in, a)

if imag(k_in)==0
    PsiMvalue = 2^(-4)*pi*a^3*(k_in*a)^2/5/gamma(5/2)^2*hypergeom([2,5/2],[5/2,7/2,4],-k_in^2*a^2);
else
    PsiMvalue = a^2*imag(conj(k_in)*spherbessJ(0,conj(k_in)*a)*spherbessJ(1,k_in*a))/(2*real(k_in)*imag(k_in));
end

% f = @(r) abs(spherbessJ(1,k_in*r)).^2.*r.^2; PsiMvalue = integral(f,0,a);

end
