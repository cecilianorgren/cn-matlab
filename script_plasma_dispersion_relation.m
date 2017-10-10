function [chi_i chi_e1 chi_e2] = chis()
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
%

chi_i = @()()





eta_i=(o-k.*vi)./k./vti;
chi_i=(2*opi.^2./(k.^2*vti.^2))*(1+eta_i*cef(eta_i,16));

%%
