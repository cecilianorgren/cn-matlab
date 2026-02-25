function out = g(llratio)
% g(l_perp/l_par): geometry factor from Tao2011 paper.

top = 1 - llratio.^2 + llratio.^2.*acos(1./llratio).*sqrt(llratio.^2-1);
bottom=(1-llratio.^2).^2;
out = real(top./bottom);