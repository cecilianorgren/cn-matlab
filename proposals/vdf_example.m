Tmax_ad = G*M*3*m/(10*2*kB*R)
Tmax_iso = G*M*1*m/(2*2*kB*R)
 

% vs_iso = sqrt(2*kB*T/m)
%
% rc = G*M./(2*vs_iso.^2)
% 
% rc/R = (2*vs/v_esc)^2
% 
% v_esc = sqrt(2*G*M/R)

%%
syms G M R T m kB
%S = syms;

v_esc = sqrt(2*G*M/R)
vs = sqrt(2*kB*T/m)
rc = G*M./(2*vs_iso.^2)


eps = (v_esc/(2*vs))^2
