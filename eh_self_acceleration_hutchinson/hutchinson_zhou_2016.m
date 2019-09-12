
nv = 100;
v1v = linspace(0,5,nv);

% Eq. 12 integrand, hole acceleration
eq12_int = @(v1v) - 2 + 3*v1v - v1v.^2;
dPoc_int = @(v1v) eq12_int(v1v);

% Eq. 18 integrand, hole growth
eq18_int = @(v1v) v1v + v1v.^3;
dPg_int = @(v1v) eq18_int(v1v);


h = setup_subplots(4,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,v1v,dPoc_int(v1v),v1v,dPoc_int(-v1v))
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'v_1/v, both velocities in frame of EH';
legend(hca,{'v_1/v>0','v_1/v<0'})


hca = h(isub); isub = isub + 1;
plot(hca,v1v,dPg_int(-v1v),v1v,dPg_int(v1v))
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'v_1/v, both velocities in frame of EH';
legend(hca,{'v_1/v>0','v_1/v<0'})


hca = h(isub); isub = isub + 1;
plot(hca,v1v,dPg_int(-v1v),v1v,1/2*dPg_int(-v1v))
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'v_1/v, both velocities in frame of EH';
legend(hca,{'v_1','v_1-> 2v_1'})

hca = h(isub); isub = isub + 1;
plot(hca,v1v,dPg_int(-v1v),v1v,1/2*dPg_int(-v1v))
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'v_1/v, both velocities in frame of EH';
legend(hca,{'v_1','v_1-> 2v_1'})

%%
% Eq. 18
eq18 = @(m,n,Udot,v1,v) m.*n.*Udot.*(- 2 + 3*v1./v - (v1./v).^3);

n = [0.080 0.020];
veh = 27000;
vAB = [0 35000];
v1 = abs(vAB-veh);
m = [1 1];
Udot = 1;
vtrap = 20000;
eq18(m,n,Udot,v1,-sqrt(v1.^2-vtrap))

%% Do integral explicitly
units = irf_units;
phimax = 2000; % V
vph = 27000*1e3; % m/s
lpar = 50*1e3; % m
m = [units.me units.me]; % kg
q = [-units.e -units.e]; % C
n = [0.080 0.020]; % cc
vd = [0 35000]*1e3; % m/s
T = [800 200]; % eV
vt = sqrt(T*units.eV*2/units.me); % m/s

% move to EH frame
vd = vd - vph; 
vph = vph - vph;

% Set up expression for momentum transfer rates
f_phi = @(x,lpar,phimax) phimax.*exp(-(x./lpar).^2/2);
f_psi = @(phi,q) phi*q; % potential energy
f_v = @(x,v1,phimax,lpar,q,m) sqrt(v1.^2-2.*f_psi(f_phi(x,lpar,phimax),q)./m); % velocity inside EH
f_f1 = @(v,vd,vt) n*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2); % distributions outside EH


%f_P_int = @(m,Udot,v1,v)

% Set up fpr plotting
nx = 100;
x = linspace(-5*lpar,5*lpar,nx);

% Plot
h = setup_subplots(4,1);
isub = 1;

if 1 %
  hca = h(isub); isub = isub + 1;
  plot(hca,x,f_phi(x,lpar,phimax))
  hca.XLabel.String = '';
  hca.YLabel.String = '';
end
if 1 %
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  plot(hca,x,f_v(x,vd(1),phimax,lpar,q(1),m(1))*1e-6,...
           x,f_v(x,vd(2),phimax,lpar,q(2),m(2))*1e-6)
  hold(hca,'on')
  plot(hca,x,f_v(x,vd(1)+vt(1),phimax,lpar,q(1),m(1))*1e-6,'--',...
           x,f_v(x,vd(1)-vt(1),phimax,lpar,q(1),m(1))*1e-6,'--')
  plot(hca,x,f_v(x,vd(2)+vt(2),phimax,lpar,q(2),m(2))*1e-6,'-.',...
           x,f_v(x,vd(2)-vt(2),phimax,lpar,q(2),m(2))*1e-6,'-.')
  hold(hca,'off')
  c_eval('hca.Children(?).Color = colors(2,:);',[1 2 5])
  c_eval('hca.Children(?).Color = colors(1,:);',[4 3 6])
  hca.XLabel.String = '';
  hca.YLabel.String = '';
  legend(hca,{sprintf('v_d = %.0f km/s',vd(1)*1e-3),sprintf('v_d = %.0f km/s',vd(2)*1e-3)})
end
if 0 %
  hca = h(isub); isub = isub + 1;
  plot(hca,[],[])
  hca.XLabel.String = '';
  hca.YLabel.String = '';
end
if 0 %
  hca = h(isub); isub = isub + 1;
  plot(hca,[],[])
  hca.XLabel.String = '';
  hca.YLabel.String = '';
end
if 0 %
  hca = h(isub); isub = isub + 1;
  plot(hca,[],[])
  hca.XLabel.String = '';
  hca.YLabel.String = '';
end
if 0 %
  hca = h(isub); isub = isub + 1;
  plot(hca,[],[])
  hca.XLabel.String = '';
  hca.YLabel.String = '';
end






















