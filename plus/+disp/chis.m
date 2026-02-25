function [chi_i chi_e1 chi_e2] = chis(o,k,R,Te1,Te2,Ti,vde1,vde2,vdi)
%% Plasma dispersion relation
% 1 + chi_i + chi_e1 + chi_e2 = 0

% needed input for chi_i:
%
% omega_pi => opi
% v_ti => vti
% k 
% omega => o
% 
% eta_i = ( omega - k*v_i)/k*v_ti
%
% Z = cef(...,16)
Te = @(R,Te1,Te2)(R*Te2+(1-R)*Te1);
chi_i = @(o,k,me,mi,Te,Ti,vdi)((2*Te/Ti/k^2)*(1 + sqrt(me*Te/mi/Ti)*o/k-vdi)*cef(sqrt(me*Te/mi/Ti)*o/k-vdi,16));
chi_e1 = @(o,k,me,mi,Te,Te1,vdi)((2*(1-R)*Te/Ti/k^2)*(1 + sqrt(me*Te/mi/Ti)*o/k-vdi)*cef(sqrt(me*Te/mi/Ti)*o/k-vdi,16));






eta_i=(o-k.*vi)./k./vti;
chi_i=(2*opi.^2./(k.^2*vti.^2))*(1+eta_i*cef(eta_i,16));

%%
