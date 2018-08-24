function [F] = get_f_maxwellian(V,n,vt,vd)  
% F = GET_F_MAXWELLIAN(V,n,vt,vd);  
units = irf_units;
  
% Set up f0
nPop = numel(n);
f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
eval(f0_str)

F = f0(V);
end