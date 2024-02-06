units = irf_units;
time = irf_time('2015-10-16T10:33:45.238Z','utc>epochtt');
times = time + [-2:1:2]*0.03*4;
times = times(1);
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
    
  c_eval('fce = fce?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
  fce = mean(fce.data); wce = fce*2*pi;
  c_eval('re = re?.tlim(pdist.time + 0.5*0.03*[-1 1]);',ic); 
  re = mean(re.data); kre = re/(2*pi);
  n = 1;
  wk0 = 0.25*wce;
  kpar0 = 0.4*kre;
  vg0 = (wk0/kpar0)*1e3; % m/s
  %vg0 = 1e6;
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

  
  
  pitch = pdist.pitchangles(b,[0:5:180]);

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

c_eval('shading(h(?),''flat'')',1:5)
hlinks = linkprop(h,{'CLim','XLim','YLim'});
colormap(pic_colors('candy4'))
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
delete(hb(1:4))

hb = hb(5);
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
