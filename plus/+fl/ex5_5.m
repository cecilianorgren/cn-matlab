% ex5_5
p_inf = 2;
R = 1;
om = 1;
rho = 1;
p_irr = @(r) p_inf-rho*om^2*R^4/2./r./r;
p_sol = @(r) p_inf+rho*om^2*(r.^2/2-R^2);

r_sol=linspace(0,R,100);
r_irr=linspace(R,4*R,100);

plot(r_irr,p_irr(r_irr),r_sol,p_sol(r_sol),R,p_sol(R),'*')
%set(gca,'ylim',[0 p_inf])
xlabel('r/R')
ylabel('p')
title('Pressure in horizontal cross-section of tornado')