%velocity_filter

vperp = 400;
vpar = 1400; % 10000 eV
vpar = 3900; % 80000 eV
vgyro = 1400; % 10000 eV

rop = 720; % 10000 eV

% lperp =  vperp * t 
% lpar = vpar * t
% 't = t' -> lperp/vperp = lpar/vpar ->  lperp = lpar*(vperp/vpar)

vfe_distance = @(lpar,vperp,vpar) lpar*vperp/vpar;

lpar_vec = linspace(0,10000,100);
plot(lpar_vec,vfe_distance(lpar_vec,vperp,vpar))
