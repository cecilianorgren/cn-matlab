% hand in 2

w=5;
l=6;
u=3;
rho=860;
eta=1e-2;
h=2e-3; % 2 mm
g=9.8;
Q = @(angle) w*(u*h/2-rho*g*h^3*sind(angle)/12/eta);

angle=0:180;
plot(angle,Q(angle),30,Q(30),'r*',angle,Q(30))

set(gca,'xlim',[0 180])
title('Flow rate')
xlabel('inclination (^o)')
ylabel('Q [m^3/s]')

%%
tau = @(angle) eta*(rho*g*h*sind(angle)/2/eta+u/h);
ang = 30;
disp(['fx = ' num2str(tau(ang)) ' N/m2'])
disp(['Fx = fx*A = ' num2str(tau(ang)*w*l) ' N'])
disp(['P = Fx*U = ' num2str(tau(ang)*w*l*u) ' Nm/s'])
%tau(ang)*w*l
%tau(ang)*w*l*u