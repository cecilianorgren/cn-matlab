function [vphi,ephi,ephikbt] = beam_solver_addons(vph,S,Te1)
% take wi and wr from beam_solver and kind maximum growth rates, constructs
% phase velocities.
%
% Assume input omega is normalized to electron plasma frequency and that k
% is normalized to 1/l_De.

%doPlot = 0;

units = irf_units;
e = units.e;
eps0 = units.eps0;
mp = units.mp;
me = units.me;
c = units.c;

vthe1 = c*sqrt(1-1./(Te1.*e./(me*c^2)+1).^2);
%vthe2 = c*sqrt(1-1./(Te2.*e./(me*c^2)+1).^2);

%vbeam = repmat(S',1,numel(R))*vthe1;
Smat=repmat(tocolumn(S),1,size(vph,2),size(vph,3),size(vph,4));

vbeam = Smat; % norm to vthe1
vphi = vbeam-vph;
ephi = vphi.^2*1e6*units.me/2/units.e/Te1;
ephi = (vphi*vthe1).^2*1e6*units.me/2/units.e/Te1;
ephikbt = ephi/Te1;
        


