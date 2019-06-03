function varargout = dispersion_solver(inp_strs,xguess,ks,powerspectra,fe1D_tint)

% load('/Users/cecilia/MATLAB/data-matlab/polarization_paper_electron_acceleration.mat')
% load('/Users/cecilia/MATLAB/data-matlab/polarization_paper_electron_acceleration_distributions_fit.mat')
% [wr,wi,k] = paper_electron_acceleration.dispersion_solver({'fast'},100,[0.1 2 100],{xvecs,yvecs,Power});
% [wr,wi,k] = paper_electron_acceleration.dispersion_solver({'slow'},100,[0.1 2 100],{xvecs,yvecs,Power});
% [wr,wi,k] = paper_electron_acceleration.dispersion_solver({'low_fast_fit'},500,[0.06 2 100],{xvecs,yvecs,Power},{ef1D,tint_disprel_fast+[0.1 0]+-0.3});
% [wr,wi,k] = paper_electron_acceleration.dispersion_solver({'fast_fit'},500,[0.06 2 100],{xvecs,yvecs,Power},{ef1D,tint_disprel_fast+[0.1 0]+-0.3});
if not(isempty(powerspectra))
  xvecs = powerspectra{1};
  yvecs = powerspectra{2};
  Power = powerspectra{3};
  doPowerspectra = 1;
else
  doPowerspectra = 0;
end
if not(isempty(fe1D_tint))
  fe1D = fe1D_tint{1};
  fe_tint = fe1D_tint{2};
  doFe = 1;
else
  doFe = 0;
end
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt.^2); % this is reduced

%% Collect input to solver
ninp_strs = numel(inp_strs);
n = [];
T = [];
m = [];
q = [];
vd = [];
vt = [];
wp = [];
Lin = [];
Ld = [];
for iinp = 1:ninp_strs
  [n_tmp,T_tmp,m_tmp,q_tmp,vd_tmp,vt_tmp,wp_tmp,Lin_tmp,Ld_tmp] = f_inp(inp_strs{iinp});
  n = [n n_tmp];
  T = [T T_tmp];
  m = [m m_tmp];
  q = [q q_tmp];
  vd = [vd vd_tmp];
  vt = [vt vt_tmp];
  wp = [wp wp_tmp];
  Lin = [Lin Lin_tmp];
  Ld = [Ld Ld_tmp];
end
nsp = numel(n); % number of species


%% Dispersion solver, one surface
% Wavenumber vector
nk = ks(3);
k_min= ks(1); k_max = ks(2);
knorm = min(Ld(1));  % length
knorm_str = sprintf('L_{d%g}',1);
kvec = linspace(k_min,k_max,nk)/knorm;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
x = xguess;

for ik = 1:nk  
  xguess = x;
  %xguess = vd(2)*kvec(ik);
  
  af = @(temp) D_streaming(temp,kvec(ik),vt,wp,vd);   
  options = optimset('GradObj','on','display','off','TolFun',1e-12);  
  [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
  fprintf('x = %g + %gi \n',real(x),imag(x))
  fval_store(ik) = FVAL; 
  wr_store(ik) = real(x);
  wi_store(ik) = imag(x);  
end

vph_store = wr_store./kvec;

rem_ind = nk;350:nk;
wr_store(rem_ind) = NaN;
wi_store(rem_ind) = NaN;
fval_store(rem_ind) = NaN;
vph_store(rem_ind) = NaN;

ikmax = find(wi_store==max(wi_store),1,'first');
kmax = kvec(ikmax);
vphmax = wr_store(ikmax)/kvec(ikmax);
wimax = wi_store(ikmax);
wrmax = wr_store(ikmax);

varargout{1} = wr_store;
varargout{2} = wi_store;
varargout{3} = kvec;

%% Plot
do_normf = 0;
wnorm = wp(1);
wnorm_str = 'w_{p1}';
h = setup_subplots(3,2,'vertical');
isub = 1;
if 1 % input distributions
  hca = h(isub); isub = isub + 1;
  vmax = max(2*vt + vd);
  vmin = min(-2*vt + vd);
  vvec = linspace(vmin,vmax,1000); 
  
  plot(hca,0,0)
  hold(hca,'on')
  % group ions and electron separately 
  fetot = 0;
  ise = find(q < 0);
  for isp = ise, fetot = fetot + f(vvec,n(isp),vt(isp),vd(isp)); end
  
  fitot = 0;
  isi = find(q > 0);
  for isp = isi, fitot = fitot + f(vvec,n(isp),vt(isp),vd(isp)); end
  
  % normalization  
  femaxtot = max(fetot);
  fimaxtot = max(fitot);
  fnorm = ones(nsp,1);
  fenormtot = 1;
  finormtot = fimaxtot/femaxtot;
  fnorm(isi) = finormtot;
  
  colors = pic_colors('matlab');
  icount = 0;
  ecount = 0;
  for isp = 1:nsp
    if q(isp) < 0, ecount = ecount + 1; displayname = sprintf('f_{e%g}',ecount); end
    if q(isp) > 0, icount = icount + 1; displayname = sprintf('f_{i%g}',icount); end
    hp_tmp = patch(hca,[vvec vvec(end) vvec(1)]*1e-6,[f(vvec,n(isp),vt(isp),vd(isp)) 0 0]/fnorm(isp),colors(isp,:),'displayname',displayname);
    hp_tmp.FaceAlpha = 0.3;
    hp_tmp.EdgeColor = 'none';    
    hp(isp) = hp_tmp;
    %f_legends{isp} = sprintf('f_%.0f/%g',isp,fnorm(isp));
  end  
 
  
  h_fe = [];
  if not(ise == 0)
    h_fe = plot(hca,vvec*1e-6,fetot/fenormtot,'-','linewidth',1,'color',[0 0 0],'displayname','f_{e,tot}');
  end
  h_fi = [];
  if not(isi == 0)
    h_fi = plot(hca,vvec*1e-6,fitot/finormtot,'-','linewidth',1,'color',[0 0 0],'displayname','f_{i,tot}');
  end
  
  % phase velocity at maximum growth rate
  h_vph = plot(hca,vphmax*1e-6*[1 1],hca.YLim,'--','linewidth',1.0,'color',0.3+[0 0 0],'displayname','v_{ph}');
  f_legends = cell(nsp,1);
 
  hleg = legend([h_fe hp(ise) h_fi hp(isi) h_vph],'Location','best');
  hleg.Box = 'off';
  
  %hl = legend(hca,{'v_{ph} @ max w_i',sprintf('f_{tot}/%g',max(ftot)),f_legends{:}},'location','best','box','off');
  %hl.Position(3) = 0.75;
  %hl.Position(2) = hl.Position(2)+0.08;
  hca.XLabel.String = 'v (10^3 km/s)';
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f_e (s/m^4)'; end
  %plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6;  
  hold(hca,'off')
end
if 1 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec*knorm,vph_store*1e-3,kmax*knorm,vphmax*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end])*knorm;
  legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = sprintf('k%s',knorm_str);
end
if 1 % solution, wr wi
  hca = h(isub); isub = isub + 1;    
  if 0
    knorm_tmp = knorm;
    knorm = 1e3;
               
    %ax = plotyy(hca,kvec*knorm,wr_store/wnorm,kvec*knorm,wi_store/wnorm);
    ax = plotyy(hca,kvec*knorm,wr_store/2/pi,kvec*knorm,5*wi_store/2/pi);  
    %ax = plotyy(hca,kvec,wr_store/wnorm,kvec*knorm,wi_store/wnorm);
    %ax(2).YLim(1) = 0;
    ax(2).YGrid = 'on';
    ax(1).XLim = kvec([1 end])*knorm;
    ax(2).XLim = kvec([1 end])*knorm;
    %ax(2).YLim = [-0.03 0.01];
    %ax(2).YTick = [-0.03:0.01:0.01];
    %ax(1).YLim = [0 1.5];
    knorm = knorm_tmp;
  else
    knorm = 1e3;
    plot(hca,kvec*knorm,wr_store/2/pi,kvec*knorm,wi_store/2/pi,'linewidth',1.5)
    hca.XLim = [0 max(kvec)*knorm];
  end
  legend(hca,'w_r','w_i', 'location','northeast')
  hca.XLabel.String = sprintf('k%s',knorm_str);
  if wnorm ~= 1
    hca.YLabel.String = sprintf('w/%s',wnorm_str);
  else
    hca.YLabel.String = '2\pi f (2\pi Hz)';
  end
  hca.YLabel.String = 'f (kHz)';
  hca.XLabel.String = 'k (km^{-1})';
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
end
if 1 % solution, wr wi  
  hca = h(isub); isub = isub + 1;     
  %%
  pcolor(hca,xvecs.kmag*1e3,yvecs.fkmag,log10(Power.Powerkmagf));
  shading(hca,'flat')
  knorm_tmp = knorm;
  knorm = 1e3;
  hold(hca,'on')
  hlines = plot(hca,kvec*knorm,wr_store/2/pi,'k-',kvec*knorm,5*wi_store/2/pi,'k--','linewidth',1.5);  
  hold(hca,'off')
  hca.XLim = kvec([1 end])*knorm;  
  hca.XLim = [0 0.25]; %km^-1
  hca.YLim(1) = 0;
  hca.YLim(2) = 1000;
  knorm = knorm_tmp;
  legend(hlines,'f_r','5f_i', 'location','northeast')
  hca.XLabel.String = 'k_{||} (km^{-1})';  
  hca.YLabel.String = 'f (Hz)';  
  
end
if doFe
  hca = h(isub); isub = isub + 1;
  fe1D.data(fe1D.data<1e-7) = NaN;
  irf_spectrogram(hca,fe1D.specrec('10^3 km/s'))
  irf_pl_mark(hca,fe_tint.epochUnix,'k')
end
if 1 % input distributions
  hca = h(isub); isub = isub + 1;
  vmax = max(2*vt + vd);
  vmin = min(-2*vt + vd);
  vvec = linspace(vmin,vmax,1000); 
  
  plot(hca,0,0)
  hold(hca,'on')
  % group ions and electron separately 
  fetot = 0;
  ise = find(q < 0);
  for isp = ise, fetot = fetot + f(vvec,n(isp),vt(isp),vd(isp)); end
  
  fitot = 0;
  isi = find(q > 0);
  for isp = isi, fitot = fitot + f(vvec,n(isp),vt(isp),vd(isp)); end
  
  % normalization  
  femaxtot = max(fetot);
  fimaxtot = max(fitot);
  fnorm = ones(nsp,1);
  fenormtot = 1;
  finormtot = fimaxtot/femaxtot;
  fnorm(isi) = finormtot;
  
  colors = pic_colors('matlab');
  icount = 0;
  ecount = 0;
  for isp = 1:nsp
    if q(isp) < 0, ecount = ecount + 1; displayname = sprintf('f_{e%g}',ecount); end
    if q(isp) > 0, icount = icount + 1; displayname = sprintf('f_{i%g}',icount); end
    hp_tmp = patch(hca,[vvec vvec(end) vvec(1)]*1e-6,[f(vvec,n(isp),vt(isp),vd(isp)) 0 0]/fnorm(isp),colors(isp,:),'displayname',displayname);
    hp_tmp.FaceAlpha = 0.3;
    hp_tmp.EdgeColor = 'none';    
    hp(isp) = hp_tmp;
    %f_legends{isp} = sprintf('f_%.0f/%g',isp,fnorm(isp));
  end  
 
  
  h_fe = [];
  if not(ise == 0)
    h_fe = plot(hca,vvec*1e-6,fetot/fenormtot,'-','linewidth',1,'color',[0 0 0],'displayname','f_{e,tot}');
  end
  h_fi = [];
  if not(isi == 0)
    h_fi = plot(hca,vvec*1e-6,fitot/finormtot,'-','linewidth',1,'color',[0 0 0],'displayname','f_{i,tot}');
  end
  
  % phase velocity at maximum growth rate
  h_vph = plot(hca,vphmax*1e-6*[1 1],hca.YLim,'--','linewidth',1.0,'color',0.3+[0 0 0],'displayname','v_{ph}');
  f_legends = cell(nsp,1);
 
  hleg = legend([h_fe hp(ise) h_fi hp(isi) h_vph],'Location','best');
  hleg.Box = 'off';
  
  %hl = legend(hca,{'v_{ph} @ max w_i',sprintf('f_{tot}/%g',max(ftot)),f_legends{:}},'location','best','box','off');
  %hl.Position(3) = 0.75;
  %hl.Position(2) = hl.Position(2)+0.08;
  hca.XLabel.String = 'v (10^3 km/s)';
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f_e (s/m^4)'; end
  %plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6;  
  hold(hca,'off')
  if doFe
    hold(hca,'on')
    plot(hca,fe1D.depend{1}(1,:)*1e-3,mean(fe1D.tlim(fe_tint).data,1),'displayname','f_{e,obs}')
    hold(hca,'off')
  end
end
%plot(kvec,wr_store,kvec,wi_store)