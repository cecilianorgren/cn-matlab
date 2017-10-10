function dispersion_relation = hot_buneman(om,k,ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi)
% Calculate the dielectric tensor.
%   disp = disp.hot_buneman(om,k,ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi);
%
Z = @(x)faddeeva(x); % Plasma dispersion function
n = @(om,k,vt,vd) (om-k*vd)./k/vt;
X = @(om,k,op,vt,vd) op^2./k.^2/vti^2*2.*(1+n(om,k,vt,vd).*Z(n(om,k,vt,vd)));

Xi = X(om,k,opi,vti,vdi);
Xe1 = X(om,k,ope1,vte1,vde1);
Xe2 = X(om,k,ope2,vte2,vde2);

dispersion_relation = 1 + Xi + Xe1 + Xe2;

disp(['Xi  = ' num2str(Xi)])
disp(['Xe1 = ' num2str(Xe1)])
disp(['Xe2 = ' num2str(Xe2)])
disp(['1 + Xi + Xe1 + Xe2 = ' num2str(dispersion_relation)])


D = @(om,k,ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi) ...
    1 + X(om,k,opi ,vti ,vdi ) + ...
        X(om,k,ope1,vte1,vde1) + ...
        X(om,k,ope2,vte2,vde2);
    
dispersion_relation = D(om,k,ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi); 
