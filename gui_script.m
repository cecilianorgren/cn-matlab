% script for running guis

% gui_changeplotparameter
% lengths in km
omega=@(k,ope,oce)((3e5)^2*k.^2*oce/ope^2);
ope=[0.01 10]*1e3; % Hz
oce=[0.01 10]*1e3; % Hz
omegaref=5e3;
k=linspace(0,1e-1,100);
gui_changeplotparameter(omega,k,ope,oce,omegaref)