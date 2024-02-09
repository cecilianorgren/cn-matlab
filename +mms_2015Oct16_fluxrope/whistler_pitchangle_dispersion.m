units = irf_units;
time = irf_time('2015-10-16T10:33:45.238Z','utc>epochtt');
times = time + [-2:1:2]*0.03*4;
times = times(5);
ic = 1;
nrows = 1;
ncols = times.length;
ip = 0;
h = gobjects([nrows,ncols]);
for irow = 1:nrows
  for icol = 1:ncols
    ip = ip + 1;
    h(ip) = subplot(nrows,ncols,ip);
  end
end



for it = 1:times.length
  time = times(it);
  c_eval('pdist = ePDist?.resample(time,''nearest'');',ic)
  c_eval('b = dmpaB?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
  c_eval('scpot = scPot?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
  
  C = @(vperp,vpar,n,wce,wk0,kpar0,vg0) vperp.^2 + n*wce./(n*wce-wk0+kpar0*vg0).*(vpar-vg0).^2;
  C = @(vperp,vpar,n,wce,wk0,kpar0,vg0) vperp.^2 + (vpar-vg0).^2;
    
  c_eval('fce = fce?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
  fce = mean(fce.data); wce = fce*2*pi;
  c_eval('re = re?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
  re = mean(re.data)*1e3; % m
  kre = (2*pi)/re;
  n = -1;
  wk0 = 0.25*wce;
  kpar0 = 0.4*kre;
  vg0 = 0.5*(wk0/kpar0); % m/s, roughly from dispersion surfaces
  vg0 = 2e6;
  Eg0 = units.me*vg0.^2/2/units.eV; % m/s
  energy = pitch.depend{1} - mean(scpot.data);
  energy(energy<0) = 0;
  energy = [0.0 energy];
  vpar = sqrt(energy*units.eV*2/units.me); % m/s
  %energy = pitch.depend{1};
  vperp = sqrt(energy*units.eV*2/units.me); % m/s
  vpar = [-vperp(end:-1:1), vperp];
  [VPERP,VPAR] = ndgrid(vperp,vpar);
  EPERP = VPERP.^2*units.me/2/units.eV;
  EPAR = sign(vpar).*VPAR.^2*units.me/2/units.eV;
  const = C(VPERP,VPAR,n,wce,wk0,kpar0,vg0);
  const_energy = const*(units.me/2/units.eV);
  const_lev = energy/(units.me/2/units.eV);

  
  
  pitch = pdist.deflux.pitchangles(b,[0:5:180]);

  hca = h(it);
  pitch.plot_pad_polar(hca,'scpot',scpot)
  clim = hca.CLim;
  if 1
    %%
    hold(hca,'on')
    %plot(hca,vpar,vperp2(const,n,wce,wk0,kpar0,vg0,vpar))
    [c,hc] = contour(hca,(EPERP),(EPAR),const,const_lev,'linewidth',0.5,'color','k');
    %clabel(c,hc)
    hold(hca,'off')
    hca.XLim = [0 800];
    hca.YLim = [-800 800];
  end
  hca.Title.String = time.utc('HH:MM:SS:mmm');
  hca.CLim = clim;
end

c_eval('shading(h(?),''flat'')',1:numel(h))
hlinks = linkprop(h,{'CLim','XLim','YLim'});
colormap(pic_colors('candy4'))
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
delete(hb(1:end-1))

hb = hb(end);
hb.Location = 'manual';
drawnow
hb.Position(1) = h(end).Position(1) + h(end).Position(3) + 0.005;

%% Theory
% Velocity-space Signatures of Resonant Energy Transfer between Whistler Waves and Electrons in the Earth's Magnetosheath
% https://iopscience.iop.org/article/10.3847/1538-4357/ad0df8
% vperp^2 + n*wce/(n*wce-wk0+kpar0*vg0)*(vpar-vg0)^2 = const.
% vpar = vg0 +- sqrt((const - vperp^2)*(n*wce-wk0+kpar0*vg0)/(n*wce))
% vperp^2 = const - n*wce/(n*wce-wk0+kpar0*vg0)*(vpar-vg0)^2
% vperp2 = @(C,n,wce,wk0,kpar0,vg0,vpar) C - n*wce/(n*wce-wk0+kpar0*vg0)*(vpar-vg0)^2;

ic = 1; 
c_eval('b = dmpaB?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
c_eval('scpot = scPot?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
c_eval('fce = fce?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
fce = mean(fce.data); wce = fce*2*pi;
c_eval('re = re?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
re = mean(re.data); kre = re/(2*pi);
n = 1;
wk0 = 0.25*wce;
kpar0 = 0.4*kre;
vg0 = wk0/kpar0;

C = @(vperp,vpar,n,wce,wk0,kpar0,vg0) vperp.^2 + n*wce./(n*wce-wk0+kpar0*vg0).*(vpar-vg0).^2;

energy = pitch.depend{1};
vperp = sqrt(energy*units.eV*2/units.me); % m/s
vpar = [-vperp(end:-1:1), vperp];
[VPERP,VPAR] = ndgrid(vperp,vpar);
const = C(VPERP,VPAR,n,wce,wk0,kpar0,vg0);

[c,h] = contour(VPERP,VPAR,const);
clabel(c,h)

%vperp2 = @(C,n,wce,wk0,kpar0,vg0,vpar) C - n*wce./(n*wce-wk0+kpar0*vg0).*(vpar-vg0).^2;
%C = @(vperp,n,wce,wk0,kpar0,vg0,vpar) vperp.^2 + n*wce./(n*wce-wk0+kpar0*vg0).*(vpar-vg0).^2;

%vperp2 = linspace(-);

%% Diffusion time
% Diffusion coefficient
units = irf_units;
dB = 0.1*1e-9; % T
fce = 555; % Hz, from local data
re = 1.0877e+03; % m, from local data
kre = (2*pi)/re; % from whamp, assuming k-normalization to be rho_e
wk0 = 0.25*wce; % from whamp
kpar0 = 0.4*kre;
vg0 = 0.5*(wk0/kpar0); % m/s, roughly from dispersion surfaces
df = 300; % eyed from spectra, but varies with time
n = -1;
vres = (wk0-wce*n)/kpar0;
Dth = (units.e/units.me)^2*dB^2/(df)*vg0/vres;

%% Diffusion time from pitchangle distribution
% dfdt = 1/sin(th)*d/dth(D*sin(th)*df/dth)

time = irf_time('2015-10-16T10:33:45.238Z','utc>epochtt');
times = time + [-2:1:2]*0.03*4;
time = times(5);
c_eval('pdist = ePDist?.resample(time,''nearest'');',ic)
c_eval('b = dmpaB?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
c_eval('scpot = scPot?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
  

%pdist = pdist.elim([0 1e4]);
dt = 5;
th_edges = 0:dt:180;
dth = diff(th_edges);
th = th_edges(1:end-1) + 0.5*dth;
pitch = pdist.pitchangles(b.norm,th_edges);
%pitch = pitch.convertto('s^3/m^6');
energy = pitch.depend{1};
energy_edges = [energy(1)-pitch.ancillary.delta_pitchangle_minus(1) energy+pitch.ancillary.delta_energy_plus];
[E,TH] = ndgrid(energy,th);

f = squeeze(pitch.data);
dfdth = f*0; 
dfdth(:,2:end-1) = (f(:,3:end)-f(:,1:end-2))/(2*dt);
dfdth(:,1) = (f(:,2)-f(:,1))/(dt);
dfdth(:,end) = (f(:,end)-f(:,end-1))/(dt);
%[~,dfdth_] = gradient(f',E,TH);
Dsinthdfdth = Dth*sind(TH).*dfdth;
%imagesc(dfdth')

%[~,dfdt__] = gradient(Dsinthdfdth,E',TH');

dfdt__ = dfdth*0; 
dfdt__(:,2:end-1) = (Dsinthdfdth(:,3:end)-Dsinthdfdth(:,1:end-2))/(2*dt);
dfdt__(:,1) = diff(Dsinthdfdth(:,1:2),1,2)/(dt);
dfdt__(:,end) = diff(Dsinthdfdth(:,end-1:end),1,2)/(dt);
dfdt = dfdt__./sind(TH);

dist_df = pitch;
dist_df.data = reshape(dfdt,[1 size(dfdt)]);
pitch_dfdt = pitch.pitchangle_diffusion(Dth);

nrows = 4;
ncols = 2;
ip = 0;
h = gobjects([nrows,ncols]);
for irow = 1:nrows
  for icol = 1:ncols
    ip = ip + 1;
    h(ip) = subplot(nrows,ncols,ip);
  end
end

isub = 1;
cmap = irf_colormap('waterfall');
if 1
  hca = h(isub); isub = isub + 1;
  pitch.deflux.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s')
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hca.XLim = [0 20];
  hca.YLim = [-10 10];
  shading(hca,'flat')
  colormap(hca,cmap)
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,energy,th,log10(f)')
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'log_{10}f';
  colormap(hca,cmap)
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,energy,th,f')
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'f';
  colormap(hca,cmap)
  %hca.CLim = prctile(abs(hca.Children.CData(:)),99)*[0 1];
  %hca.CLim = max(abs(hca.Children.CData(:)))*[-0 1];
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,energy,th,dfdth')
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'df/dth';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.Children.CData(:)))*[-1 1];
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,energy,th,Dsinthdfdth')
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'D*sin(th)*df/dth';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.Children.CData(:)))*[-1 1];
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,energy,th,dfdt')
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'df/dt';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.Children.CData(:)))*[-1 1];
end
if 1
  hca = h(isub); isub = isub + 1;
  [~,~,hh] = dist_df.plot_pad_polar(hca,'scpot',scpot,'nolog10','flim',[-Inf Inf],'10^3 km/s');
  hb = hh.Colorbar;
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hca.XLim = [0 20];
  hca.YLim = [-10 10];
  shading(hca,'flat')
  hb.YLabel.String = 'df/dt';
  hca.CLim = prctile(abs(hca.Children.CData(:)),95.)*[-1 1];
  %hca.CLim = max(abs(hca.Children.CData(:)))*[-1 1];
  colormap(hca,pic_colors('blue_red'))
  hold(hca,'on')
  plot(hca,hca.XLim,vres*1e-6*[1 1],'k--')
  hold(hca,'off')
end
if 0
  hca = h(isub); isub = isub + 1;
  toplot = (f./dfdt)';
  %toplot(isnan(toplot)) = 0;
  pcolor(hca,energy,th,toplot)
  hca.XScale = 'log';
  shading(hca,'flat')
  hca.XLabel.String = 'E';
  hca.YLabel.String = '\theta';
  hb = colorbar(hca);
  hb.YLabel.String = 'f/(df/dt)';
  colormap(hca,pic_colors('blue_red'))
  %hca.CLim = max(abs(hca.Children.CData(:)))*[-1 1];
  hca.CLim = prctile(abs(hca.Children.CData(:)),95.)*[-1 1];
  hold(hca,'on')
  contour(hca,energy,th,smooth2((f),3)','k')
  hold(hca,'off')
end

if 1
  hca = h(isub); isub = isub + 1;
  pitch_dfdt.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s','flim',[-inf inf],'nolog10')
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hca.XLim = [0 20];
  hca.YLim = [-10 10];
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = prctile(abs(hca.Children.CData(:)),95.)*[-1 1];
end




%hlinks = linkprop(h([1 2]),{'CLim'});
%h(1).CLim = [-31 -25];

hlinks = linkprop(h([1 7]),{'XLim','YLim'});
%h(1).CLim = [-31 -25];

hlinks = linkprop(h([6 7]),{'CLim'});

%% Diffusion time from pitchangle distribution
% dfdt = 1/sin(th)*d/dth(D*sin(th)*df/dth)

time = irf_time('2015-10-16T10:33:45.238Z','utc>epochtt');
%times = time + [-2:1:2]*0.03*4+0.03*4*2;
times = time + [-2:1:2]*0.03*4;
%time = times(5);


%times = times(4);
nrows = 4;
ncols = times.length;
%nrows = times.length;
%ncols = 4;

ip = 0;
h = gobjects([nrows,ncols]);
for irow = 1:nrows
  for icol = 1:ncols
    ip = ip + 1;
    h(irow,icol) = subplot(nrows,ncols,ip);
  end
end
h = h;
dt = 7.5;
th_edges = 0:dt:180;

for it = 1:times.length
  time = times(it);
  c_eval('pdist = ePDist?.resample(time,''nearest'');',ic)
  c_eval('b = dmpaB?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
  c_eval('scpot = scPot?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic)
  z = mean(b.norm.data,1);
  x = cross(z,cross([1 0 0],x)); x = x/norm(x);
  y = cross(z,x); y = y/norm(y);
  vg = -15000:500:15000;
  
  pitch = pdist.pitchangles(b.norm,th_edges);
  pitch_dfdt = pitch.pitchangle_diffusion(Dth);
  pdist_red = pdist.reduce('2D',x,z,'vg',vg);
  
  
  vz = -3000; % km/s
  vshift = z*vz;
  
  new_E = logspace(-1,log10(40e3),50);
  new_az = 0:11.25:360;
  new_pol = 0:11.25:180;
  
  %new_E = [pdist.depend{1}(1)-pdist.ancillary.delta_energy_minus(1) pdist.depend{1}-+pdist.ancillary.delta_energy_minus];
  %new_az = pdist.depend{2};
  %new_pol = pdist.depend{3};
  pdist_shift = pdist.shift(vshift,1000,[1 0 0; 0 1 0; 0 0 1],'','new_grid',new_E,new_pol,new_az);
  %pdist_shift = pdist.shift(vshift,1000,[1 0 0; 0 1 0; 0 0 1],'');
  pdist_shift = pdist_shift.convertto('s^3/cm^6');
  
  
  pitch_vres = pdist_shift.pitchangles(b.norm,th_edges);
  pitch_vres_dfdt = pitch_vres.pitchangle_diffusion(Dth);
  pdist_vres_red = pdist_shift.reduce('2D',x,z,'vg',vg);

  isub = 1;
  cmap = irf_colormap('waterfall');
  if 0 % def
    hca = h(isub,it); isub = isub + 1;
    pitch.deflux.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s','nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [0 20];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,cmap)
    hca.Title.String = pdist.time.utc('yyyy-mm-ddTHH:MM:SS.mmm');
  end
  if 0 % def shift
    hca = h(isub,it); isub = isub + 1;
    pitch_vres.deflux.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s','nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [0 20];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,cmap)
    hca.Title.String = pdist.time.utc('yyyy-mm-ddTHH:MM:SS.mmm');
  end
  if 1 % f
    hca = h(isub,it); isub = isub + 1;
    pitch.deflux.plot_pad_polar(hca,'scpot',scpot.resample(pitch),'10^3 km/s','nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [0 10];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,cmap)
    
    hold(hca,'on')
    plot(hca,hca.XLim,(vres*1e-6-0)*[1 1],'k-')
    plot(hca,hca.XLim,(vz*1e-3-0)*[1 1],'k--')
    hold(hca,'off')
  end
  if 1 % f shift
    hca = h(isub,it); isub = isub + 1;
    pitch_vres.plot_pad_polar(hca,'10^3 km/s','nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [-0 10];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,cmap)
    hold(hca,'on')
    plot(hca,hca.XLim,(vres*1e-6-vz*1e-3)*[1 1],'k-')
    plot(hca,hca.XLim,(vz*1e-3-vz*1e-3)*[1 1],'k--')
    hold(hca,'off')
  end
  if 0 % f reduced
    hca = h(isub,it); isub = isub + 1;
    pdist_red.plot_plane(hca,'log10',1);
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [-0 15];
    hca.YLim = [-15 15];
    shading(hca,'flat')
    colormap(hca,cmap)
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
  end
  if 0 % f reuced shifted
    hca = h(isub,it); isub = isub + 1;
    pdist_vres_red.plot_plane(hca,'log10',1);
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [-0 15];
    hca.YLim = [-15 15];
    shading(hca,'flat')
    colormap(hca,cmap)
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
  end
  if 1 % df/dt
    hca = h(isub,it); isub = isub + 1;
    pitch_dfdt.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s','flim',[-inf inf],'nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [-0 10];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = prctile(abs(hca.Children.CData(:)),95.)*[-1 1];
    hold(hca,'on')
    plot(hca,hca.XLim,(vres*1e-6-0)*[1 1],'k-')
    plot(hca,hca.XLim,(vz*1e-3-0)*[1 1],'k--')
    hold(hca,'off')
  end
  if 1 % df/dt
    hca = h(isub,it); isub = isub + 1;
    pitch_vres_dfdt.plot_pad_polar(hca,'scpot',scpot,'10^3 km/s','flim',[-inf inf],'nolog10')
    %hca.XScale = 'log';
    %hca.YScale = 'log';
    hca.XLim = [-0 10];
    hca.YLim = [-10 10];
    shading(hca,'flat')
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = prctile(abs(hca.Children.CData(:)),95.)*[-1 1];
    hold(hca,'on')
    plot(hca,hca.XLim,(vres*1e-6-vz*1e-3)*[1 1],'k-')
    plot(hca,hca.XLim,(vz*1e-3-vz*1e-3)*[1 1],'k--')
    hold(hca,'off')
  end
end



c_eval('hlinks? = linkprop(h(?,:),{''CLim''});',1:nrows)
%hlinks_all = linkprop(h(:),{'XLim','YLim'});
%h(1).XLim = [0 10];
%hlinks_all = linkprop(h(1:2),{'XLim','YLim','CLim'});
hlinks_1a = linkprop(h(1,:),{'XLim','YLim','CLim'});
hlinks_1b = linkprop(h(2,:),{'XLim','YLim','CLim'});
%h(1).CLim = [0 7]*1e-26;
hlinks_2 = linkprop(h(3:4),{'XLim','YLim','CLim'});

% for ip = 1:numel(h)
%   hca = h(ip);
%   hold(hca,'on')
%   plot(hca,hca.XLim,vres*1e-6*[1 1],'k--')
%   plot(hca,hca.XLim,-vshift(3)*1e-3*[1 1],'k--')
%   hold(hca,'off')
% end
drawnow
%compact_panels(h,0.01,0.005)
%hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
%delete(hb(1:end-4))
%for ip = (ncols):numel(h)
%  hca = h(ip);
%  hca.YLabel = [];
%  hca.YTick = [];
%end