elim = [300 Inf];
M1 = iPDist3.elim(elim).moments;
M2 = mms.psd_moments(iPDist3,scPot3,'energyrange',elim);
MP = iPDist3.elim(elim).macroparticles('nbin',1);
M5 = pdist_cleaned.elim(elim).moments;
N6 = iPDist3.elim(elim).n;

if 1
  %% Wrong by seemingly exactly a factor 4*pi
  n = zeros(iPDist3.length,1);
  vel = zeros(iPDist3.length,3);
  
  for it = 1:iPDist3.length
    dv = MP(it).dv;
    df = MP(it).df;
    djx = MP(it).vx.*dv.*df;
    djy = MP(it).vy.*dv.*df;
    djz = MP(it).vz.*dv.*df;

    n(it) = sum(dv.*df);
    vel(it,1) = sum(djx)./n(it);
    vel(it,2) = sum(djy)./n(it);
    vel(it,3) = sum(djz)./n(it);
  end
  M3.n = irf.ts_scalar(iPDist3.time,n);
  M3.v = irf.ts_vec_xyz(iPDist3.time,vel);
end

% M4.n = irf.ts_scalar(iPDist3.time,nn);
irf_plot({M1.n,M2.n_psd,M3.n,M5.n,N6},'comp');

%%
PD = iPDist3;
elim = [300 Inf];
M1 = PD.elim(elim).moments;
M2 = mms.psd_moments(PD,scPot3*0,'energyrange',elim);
N3 = PD.elim(elim).n;
V3 = PD.elim(elim).vel;

figure; irf_plot({M1.n,M2.n_psd,N3},'comp');
figure; irf_plot({M1.V,M2.V_psd,V3},'comp');
figure; irf_plot({PD.elim([0 Inf]).vel, PD.elim([500 Inf]).vel, PD.elim([1000 Inf]).vel},'comp')
figure; irf_plot({PD.elim([0 Inf]).vel.x, PD.elim([500 Inf]).vel.x, PD.elim([1000 Inf]).vel.x},'comp')
%% Plot

h = irf_plot(1);

hca = irf_panel('density');
irf_plot(hca,{M1.n,M2.n_psd,M3.n});
hca.YLabel.String = 'n (cc)';


%% Check why PDist.macroparticles  give the wrong density

%%%original coordinate system original reference frame
u = irf_units;
PD_FA = iPDist3(1);

En = PD_FA.depend{1};
Vn = sqrt(2*En*u.e/u.me);
thn = PD_FA.depend{3};
phn = PD_FA.depend{2};

E_minus = (PD_FA.depend{1} - PD_FA.ancillary.delta_energy_minus);
E_plus  = (PD_FA.depend{1} + PD_FA.ancillary.delta_energy_plus);
v_minus = sqrt(2*units.e*E_minus/units.me); % m/s
v_plus  = sqrt(2*units.e*E_plus/units.me); % m/s
V_edges = [v_minus v_plus(end)];
phi_edges = 0:11.25:360;
theta_edges = 0:11.25:180;
      

Vedges = V_edges*1000;
thedges = theta_edges;
phedges = phi_edges;

%%calculate differential elements


dv = diff(Vedges);
dth = diff(thedges)*pi/180;
dph = diff(phedges)*pi/180;


%%calculate d3v for F_3D
[VV,PH,TH] = ndgrid(Vn,phn,thn);
[dV,dPH,dTH] = ndgrid(dv,dph,dth);
d3v_1 = VV.^2.*sind(TH).*dV.*dTH.*dPH;

d3v_2 = PD_FA.d3v.data;