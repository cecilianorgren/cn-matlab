% Formation of powerlaw distribution by Fermi acceleration
% dw = v^2w
% dw/w = dln(w) = v^2    (dN?)
% w = exp(v^2N)
% Particle moving faster will reflect more times... And it will not be a
% linear increase, since the velocity is increasing...
% w(t) = exp(v^2t/tau)
v = 1;
N = 0:1:100;
loglog(N,exp(v^2*N).*exp())