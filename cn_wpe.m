function [out] = cn_wpe(in,flag)
% CN_WPE
% Calculates the electron plasma frequency in kHz (flag=1) 
% or the density in cc (flag=-1).

irf_units;
switch flag
    case 1
        n=in*1e6; % per cubic meter
        Wpe = sqrt(n*Units.e^2/Units.me/Units.eps0); % rad/s
        Fpe = Wpe/2/pi; % Hz
        out=Fpe*1e-3; % kHz
    case -1
        Fpe=in;
        Wpe=Fpe*2*pi;
        n=Wpe.^2*Units.me*Units.eps0/(Units.e^2);
        out=n*1e-6; % cc
end
