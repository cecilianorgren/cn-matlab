function [f0,params] = get_f0(inp);

switch inp 
  case 1
    ntot = 0.04*1e6;
    R = 0.60; 
    n1 = ntot*R; 
    n2 = ntot*(1-R);
    T1 = 150;
    T2 = 2000;
    vd1 = -9000*1e3;
    vd2 = 4000*1e3; 
  case 2
    ntot = 0.04*1e6;
    R = 0.60; 
    n1 = ntot*R; 
    n2 = ntot*(1-R);
    T1 = 50;
    T2 = 1000;
    vd1 = -18000*1e3;
    vd2 = 4000*1e3;     
end

n = [n1 n2];
vd = [vd1 vd2];
T = [T1 T2];

units =irf_units;
vt = sqrt(2*units.e*T./units.me); % m/s

nPop = 2;
if 0
f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
else
f0_str = ['f0 = @(v,n,vd,vt) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
end
eval(f0_str)

params.n = [n1 n2];
params.vd = [vd1 vd2];
params.T = [T1 T2];
params.vt = vt;