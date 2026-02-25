% Buneman1959, he assumes Te = Ti
units = irf_units;
u_min = @(Te,Ti) 0.926*(sqrt(2*units.kB*Te*11604/units.me)+sqrt(2*units.kB*Ti*11604/units.mp));

Te = 0.5:0.1:1.5; Te = linspace(0.5,2,10); nTe = numel(Te);
Ti = 0.5:0.1:1.5; Ti = linspace(0.5,2,20); nTi = numel(Ti);

u_crit = zeros(nTe,nTi);

for ie=1:nTe
    for ip=1:nTi
        u_crit(ie,ip) = u_min(Te(ie),Ti(ip));
    end
end

surf(Ti,Te,u_crit/u_min(1,1))
xlabel('Ti')
ylabel('Te')
zlabel('u_{crit}/u_{norm}')

%% Buneman1959, he assumes Te = Ti
units = irf_units;
u_thres = @(TT) 0.926*(1+sqrt(1./TT*units.me/units.mp));

% TT = Ti/Te
TT = linspace(0.1,10,100);

plot(u_thres(TT),TT)
xlabel('u/v_{te}')
ylabel('T_e/T_i')


%% Gary: Theory of space plasma microinstabilities p.35
units = irf_units;
u_min = @(Te,Ti) 1+(Te/Ti)^(3/2)*sqrt(units.mp/units.me)*exp(-1.5-Te/2/Ti);

Te = 0.5:0.1:5.5; Te = linspace(0.1,10,130); nTe = numel(Te);
Ti = 0.5:0.1:1.5; Ti = linspace(0.1,2,125); nTi = numel(Ti);

u_crit = zeros(nTe,nTi);

for ie=1:nTe
    for ip=1:nTi
        u_crit(ie,ip) = u_min(Te(ie),Ti(ip));
        c_s(ie,ip) = sqrt((3*Ti(ip)+Te(ie))*units.kB*11604/units.mp);
        v_te(ie,ip) = sqrt(2*Te(ie)*units.kB*11604/units.me);
    end
end

%surf(Ti,Te,u_crit./c_s)
surf(Ti,Te,u_crit); hold on; % .*c_s./v_te*sqrt(units.mp/units.me)
xlabel('Ti')
ylabel('Te')
zlabel('u_{crit}/v_{ph}')
ch = colorbar('peer',gca);
ylabel(ch,'u_{crit}/v_{ph}')

x = Ti;
y = 3*Ti;
plot(x,y,'*')

plot(x,x,'k','linewidth',3)



hold off

%% Deabye length ratios

DD = @(Tb,Tbg,R) sqrt(Tbg/Tb)*sqrt((1-R)./R);
Tb = 60;
Tbg = 2000;
R=linspace(0.01,1,100);
plot(R,DD(Tb,Tbg,R),R,DD(Tb,1600,R))
xlabel('R')
ylabel('d_{b}/d_{bg}')






