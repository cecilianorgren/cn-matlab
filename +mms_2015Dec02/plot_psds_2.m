%% Example plot 1
%mms_2015Oct16.plot_psd_projection
% toplot = 4; % knee in gyrotropic perp, 6 pitchangles
% toplot = 3; % 
% toplot = 6; % mva: 3pa, 3 projections

%% Example plot 2
% Define time interval
tint = irf.tint('2015-12-02T01:14:56.6Z',0.04);

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));

% Initialize figure
for ii = 1:4; h(ii) = subplot(2,2,ii); end
isub = 1;

energies = [43 178];
c_eval('dist = desDist?',ic)
% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',energies(1),'vectors',{hatB0,'B'},'flat');

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',energies(2),'vectors',{hatB0,'B'},'flat');

% Plot projection onto a plane perpendicular to B
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; hatB0;0 1 0; ],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B'});
colormap(hca,'jet')
hca.CLim = [-1 4.5];

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',energies,''ylim'',[1e-2 1e6])',ic)

if 0
  v = linspace(100,15000,100);
  units = irf_units;
  E = units.me*v.^2*1e6/2/units.e;
  f = cn.maxwellian(v,100,5,0,0,1);
  hold(hca,'on')
  loglog(hca,E,f*1e30,'--')
  
elseif 0
  hold(hca,'on')
  v = logspace(3,5,100);
  xEeV =units.me*(v*1e3).^2/2/units.e;

  nper = 25; Tper = 35; fper = cn.maxwellian(v,Tper,nper,0,'e',1);
  npar = 20; Tpar = 35; fpar = cn.maxwellian(v,Tpar,npar,0,'e',1);

  maxlines = loglog(hca,xEeV,fper*1e18,'--');%,xEeV,fpar*1e18,'--');
  %legend(maxlines,{['n = ' num2str(nper) ' cm^{-3}  T = ' num2str(Tper) ' eV'],...
  %                 ['n = ' num2str(npar) ' cm^{-3}  T = ' num2str(Tpar) ' eV']},'location','northoutside')
  legend(maxlines,['n = ' num2str(nper) ' cm^{-3}  T = ' num2str(Tper) ' eV'],'location','northoutside')
  hca.YLim = [1e-2 1e6];
  %hca.XLim = [1e-2 1e6];
  hold(hca,'off')
end
%% Other time interval
tint = irf.tint('2015-10-16T10:33:30.25Z',0.04); 

if 0 % get mva direction for projection plot
  %%
  % Automatic results made from gui above
  ic = 1:4;
  tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
  c_eval('[out?,l?,mva?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
  mva = mva3;
  % rotate around x-axis 
  turn_angle = 0; % degrees
  tmpv = mva;
  mva(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
  mva(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
  % v is 3x3 matrix:
  %v1: [0.0906 0.0896 0.9918] - first row
  %v2: [0.3510 -0.9349 0.0524] - second row
  %v3: [0.9320 0.3434 -0.1161] - third row   
end

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));
energy1 = 50;
energy2 = 180;

% Initialize figure
for ii = 1:4; h(ii) = subplot(2,2,ii); end
isub = 1;

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',energy1,'vectors',{hatB0,'B'});

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',energy2,'vectors',{hatB0,'B'});

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energylevel',[10 11 12],'vectors',{hatB0,'B'},'flat');

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',energy2,'vectors',{hatB0,'B'},'flat');

%%
% Plot projection onto a plane perpendicular to B
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,desDist1,'tint',tint,'xyz',mva,'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B'},'vlabel',{'v_L','v_M','v_N'});

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1,'energies',[energy1 energy2])


%% Example plot 3
% Define time interval
tint = irf.tint('2015-10-16T10:33:45.75Z',0.04);

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));


% Initialize figure
for ii = 1:4; h(ii) = subplot(2,2,ii); end
isub = 1;

c_eval('dist = desDist?;',ic)
% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',150,'vectors',{hatB0,'B'});

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',150,'vectors',{hatB0,'B';[0 0 1],'Z'},'flat');

% Plot projection onto a plane perpendicular to B
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; hatB0; cross(hatB0,cross(hatB0,[0 1 0]))],'elevationlim',20,'vlim',25000,'vectors',{hatB0,'B'});

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint,''scPot'',P?,''energies'',[70 180])',ic)

%% Compare fpi_plot_proj (andreas function) with mms.plot_projection

tint = irf.tint('2015-10-16T10:33:45.00Z',0.2)+0*(-6);
vlim = 15*1e3;
elevlim = 90;
correctBin = 0;


c_eval('dist = desDist?;',ic)

c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatVi0 = double(irf_norm(Vi0));
c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatVe0 = double(irf_norm(Ve0));

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));

% Initialize figure
for ii = 1:8; h(ii) = subplot(2,4,ii); end

isub = 1;

% Plot projection onto a plane
hca = h(isub);  isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; 0 0 1; 0 1 0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 1 0; 0 0 1],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; cross(hatB0,[0 1 0])],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; hatB0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);
%%
if 0
  %%
%hca = h(isub); isub = isub + 1; 
%mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[1 0 0; 0 0 1; -1.5 -1 0],'elevationlim',elevlim,'vlim',vlim,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'usebincorrection',correctBin);

%mms.plot_skymap(hca,desDist1,'tint',tint(1),'energy',10,'vectors',{hatB0,'B'},'flat','log');

hca = h(isub); isub = isub + 1; 
mms.plot_skymap(hca,dist,'tint',tint(1),'energy',1000,'vectors',{hatB0,'B';hatVi0,'V_i';hatVe0,'V_e'},'log');

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,vi?brst.tlim(tint+5*[-1 1]));',ic)
irf_pl_mark(hca,tint.epochUnix','green')
irf_zoom(hca,'x',tint+5*[-1 1])
irf_legend(hca,{'v_x','v_y','v_z'},[0.02 0.98])

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,ve?brst.tlim(tint+5*[-1 1]));',ic)
irf_pl_mark(hca,tint.epochUnix','green')
irf_zoom(hca,'x',tint+5*[-1 1])
irf_legend(hca,{'v_x','v_y','v_z'},[0.02 0.98])
end

% Andreas function
hca = h(isub); isub = isub + 1;
fpi_plot_proj(hca,dist,tint,'xy')
hca.YLim = vlim*[-1 1]; hca.XLim = hca.YLim;

hca = h(isub); isub = isub + 1;
fpi_plot_proj(hca,dist,tint(1),'xy')
hca.YLim = vlim*[-1 1]; hca.XLim = hca.YLim;
 
hca = h(isub); isub = isub + 1;
fpi_plot_proj(hca,dist,tint(1),'xz')
hca.YLim = vlim*[-1 1]; hca.XLim = hca.YLim;

%hca = h(isub); isub = isub + 1;
%mms.plot_skymap(hca,desDist1,'tint',tint(1),'energy',20,'vectors',{hatB0,'B'},'flat','log');
%%
for ii = [1 2 3 4 5 6 7];
  h(ii).CLim = [3 9];
end



%%

