%% Maxwellian
v=linspace(-1e8,1e8,1e5); vd=5e6; vt=3e7; 
maxwellian=@(v,vd,vt) exp(-((v-vd)/vt).^2)./vt/sqrt(pi);
%plot(v,maxwellian(v,vd,vt))

% To what velocity does ephi correspond
irf_units;
phi=150; % V
epot=Units.e*phi;
v_phi=sqrt(2*epot/Units.me);

% Kinetic energy
ekin=@(v)(0.5*Units.me*v.^2);
energy=@(v,vd,vt)(0.5*Units.me*(v).^2.*maxwellian(v,vd,vt));

% Integrate distribution to see how large part has kinetic energy lower
% than potential energy in wave frame
fit_density=fit(tocolumn(v),tocolumn(maxwellian(v,vd,vt)),'linearinterp');

int_density=integrate(fit_density,v,v(1));
sum_density=int_density(end)-int_density(1);

sum_density_vphi=int_density(find(v<v_phi,1,'last'))-int_density(find(v>(-v_phi),1,'first'));
fraction_trapped=sum_density_vphi/sum_density

fit_energy=fit(tocolumn(v),tocolumn(energy(v,vd,vt)),'linearinterp');
int_energy=integrate(fit_energy,v,v(1));

int_energy=integrate(fit_energy,v,v(1));
sum_energy=int_energy(2)-int_energy(1);

% Look at distributions of density and energy
maxf=max(maxwellian(v,vd,vt));
maxe=max(energy(v,vd,vt));
subplot(4,1,1);plot(v,maxwellian(v,vd,vt),[-v_phi -v_phi],[0 maxf],[v_phi v_phi],[0 maxf])
legend('f','-vphi','+vphi');ylabel('f');xlabel('v')
subplot(4,1,2);plot(v,energy(v,vd,vt),[-v_phi -v_phi],[0 maxe],[v_phi v_phi],[0 maxe],[v(1) v(end)],[epot epot])
legend('E','-vphi','+vphi');ylabel('E');xlabel('v')

% Look at distributions of integrated density and energy
maxintf=max(int_density);
maxinte=max(int_energy);
subplot(4,1,3);plot(v,int_density,[-v_phi -v_phi],[0 maxintf],[v_phi v_phi],[0 maxintf])
legend('int(f)','-vphi','+vphi');ylabel('integrated f');xlabel('v')
subplot(4,1,4);plot(v,int_energy,[-v_phi -v_phi],[0 maxinte],[v_phi v_phi],[0 maxinte])
legend('int(E)','-vphi','+vphi');ylabel('integrated E');xlabel('v')
%% Look at Maxwellian with mv^2/2 on x-axis
maxf=max(maxwellian(v,vd,vt));
plot(ekin(v),maxwellian(v,vd,vt),-[epot epot],[0 maxf],[epot epot],[0 maxf])
legend('f','-vphi','+vphi')
ylabel('f');xlabel('0.5*mev^2')

%%

%%
plot(v,int_density,v,int_density_lim)
%% To what velocity does ephi correspond
irf_units;
phi=300; % V
epot=Units.e*phi;
v_phi=sqrt(2*Units.e*phi/Units.me);

energy=@(v,vd,vt)(0.5*Units.me*v.^2.*maxwellian(v,vd,vt))
fit_energy=fit(tocolumn(v),tocolumn(energy(v,vd,vt)),'linearinterp');

int_energy=integrate(fo2,v,v(1));
sum_energy=int2(2)-int2(1);

int_energy_limited=integrate(fit_energy,-v_phi,v_phi);

%%
%%
%%
