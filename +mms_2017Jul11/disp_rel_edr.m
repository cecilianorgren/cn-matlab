% disp_rel_edr
% solver_streaming/solver_streaming.m
% Multi-species streaming instability
% Solver for 1D unmagnetized plasma
% If solution is not as expected, change initial guess: x = 0; 
% i-e instability typically low ~wpi, 
% e-e bump-on-tail instability typically higher ~vd*kvec(ik);
% e-e two-stream instablity typically vbulk*k(1)
% If solution seems ok but out of k range change knorm or k_min, k_max.

% Physical constants
units = irf_units;
qe = units.e; % 1.602e-19;
me = units.me; % 9.109e-31;
mi = units.mp; % 1.673e-27;
mime = mi/me; % 1836;
eps0 = units.eps0; % 8.854e-12;
kB = units.kB; % 1.381e-23;
c = units.c; % 299792458;
%% Out-of-plane beam direction, using a 1D beam (i.e. "parallel")
% Time for distribution fit
time = irf_time('2017-07-11T22:34:03.000Z','utc>epochtt');

% Properties for the different plasma species
B = 10e-9; % not used
n = [0.03 0.03]*1e6;
T = [300 8000]; T_K = T*units.eV/kB; % use parallel temperature
m = [me mi];
q = [-1 1]*qe; 
vd = [18000e3 0]; % m/s

nsp = numel(n); % number of species

% Physical parameters
vt = sqrt(2*qe*T./m); % m/s
wp = sqrt(n.*q.^2./(m*eps0)); % Hz
wc = q*B./m; % not used, although it can affect stability of phase space holes
ro = vt./wc;
Lin = c./wp;
Ld = vt./wp/sqrt(2);

ntot = sum(n); wptot = sqrt(ntot.*q(1).^2./(m(1)*eps0)); % Hz
toPlot = 1:2;

disp('----------------------------------------------')
fprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n))
fprintf('n = ['); fprintf(' %g',n*1e-6); fprintf('] cc\n');
fprintf('m = ['); fprintf(' %g',m/me); fprintf('] me\n');
fprintf('q = ['); fprintf(' %g',q/qe); fprintf('] e\n');
fprintf('T = ['); fprintf(' %g',T); fprintf('] eV\n');
fprintf('vt = ['); fprintf(' %.0f',vt*1e-3); fprintf('] km/s\n');
fprintf('vd = ['); fprintf(' %.0f',vd*1e-3); fprintf('] km/s\n');
fprintf('wp = ['); fprintf(' %g',wp); fprintf('] Hz\n');
fprintf('wc = ['); fprintf(' %g',wc); fprintf('] Hz\n');
fprintf('ro = ['); fprintf(' %g',ro*1e-3); fprintf('] km\n');
fprintf('Lin = ['); fprintf(' %g',Lin*1e-3); fprintf('] km\n');
fprintf('Ld = ['); fprintf(' %g',Ld*1e-3); fprintf('] km\n');
disp('----------------------------------------------')


% Plot distributions
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(3/2)*exp(-(v-vd).^2/vt.^2);
if 0 % also done below after solver
  figure(52)
  f_legends = cell(nsp,1);
  ud = get(gcf);
  delete(ud.Children)
  if 1 % plot everything normalized int eh same panel
    for isp = 1:nsp
      hca = subplot(1,1,1);
      vmax = max(vt) + max(vd);
      vvec = linspace(-2*vmax,2*vmax,1000);
      hold(hca,'on')
      plot(hca,vvec*1e-6,f(vvec,n(isp),vt(isp),vd(isp))/max(f(vvec,n(isp),vt(isp),vd(isp))))
      hold(hca,'off')
      f_legends{isp} = sprintf('f_%.0f/%g',isp,max(f(vvec,n(isp),vt(isp),vd(isp))));
    end
    legend(f_legends{:});
    hca.XLabel.String = 'v (10^3 km/s)';
  else % plot everything in different subplots
    for isp = 1:nsp
      hca = subplot(nsp,1,isp);
      vvec = linspace(-2*max(vt),2*max(vt),1000);
      plot(hca,vvec*1e-3,f(vvec,n(isp),vt(isp),vd(isp)));
      f_legends{isp} = sprintf('f_%.0f',isp);
      legend(f_legends{isp});
      hca.XLabel.String = 'v (10^3 km/s)';
      hca.YLabel.String = 'f (s^3/m^6)';
    end  
  end
end
%h = subplot(1,1,1);
if 0 % Measured electron distribution and WHAMP fit
  hca = h(isub); isub = isub +1;  
  whamp.plot_f(hca,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XLim = [1e1 5e3];
  hca.YLim = [1e-2  1e6];
  hca.XTick = [1e-1 1e0 1e1 1e2 1e3];

  % add real distribution  
  tInd = find(abs(ePitch1.time-time)==min(abs(ePitch1.time-time)));
  timeUTC = time.utc;
  hold(hca,'on')
  unitscale = 1e30; % cm^-6 ->km^-6
  mms_pa = plot(hca,ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,1)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,7)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,13)*unitscale,'+');
  mms_pa(1).Color = hca.ColorOrder(1,:);
  mms_pa(2).Color = hca.ColorOrder(2,:);
  mms_pa(3).Color = hca.ColorOrder(3,:);
  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLim = [1e-32 2e-25]*unitscale;
  hca.XLim = [1e1 3e3];  
  hca.YLabel.String = 'f_e (s^3 km^{-6})';
  hca.XTick = [1e1 1e2 1e3 1e4];
  %hleg = irf_legend(h,{'0';'90';'180'},[0.98 0.98]);
  hca.Title.String = ['Electron distribution'];
  %hold(h,'on')
  %mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])
  %hold(h,'off')
  irf_legend(hca,{{'-  Fit';'+ Measured';['   ' timeUTC(1:10)];['   ' timeUTC(12:23)]}},[0.05 0.4],'color',[0 0 0])
end


% Dispersion solver, one curve
% Wavenumber vector
nk = 700;
k_min= 0.1; k_max = 0.7;
knorm = min(Ld(1));  % length
knorm_str = sprintf('L_{d%g}',1);
kvec = linspace(k_min,k_max,nk)/knorm;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
x = 000;
% Step through k vector
for ik = 1:nk  
  xguess = x;
  %xguess = vd(2)*kvec(ik);
  
  af = @(temp) D_streaming(temp,kvec(ik),vt,wp,vd);   
  options = optimset('GradObj','on','display','off','TolFun',1e-12);  
  [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
  %fprintf('x = %g + %gi \n',real(x),imag(x));
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

%% plot solution
figure(81)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [10 10 700 1100];

wnorm = wp(1); wnorm_str = sprintf('w_{p%g}',1);

do_normf = 1;

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;



if 1 % Measured electron distribution and WHAMP fit
  hca = h(isub); isub = isub +1;  
  d = [1 1];
  a1 = [1 1];
  a2 = [0 0];
  toPlot = 1;
  whamp.plot_f(hca,n(toPlot)*1e-6,m(toPlot)/me,T(toPlot)*1e-3,vd(toPlot)./vt(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XLim = [1e1 5e3];
  hca.YLim = [1e-2  1e6];
  hca.XTick = [1e-1 1e0 1e1 1e2 1e3];

  % add real distribution  
  tInd = find(abs(ePitch1.time-time)==min(abs(ePitch1.time-time)));
  timeUTC = time.utc;
  hold(hca,'on')
  unitscale = 1e30; % cm^-6 ->km^-6
  mms_pa = plot(hca,ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,1)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,7)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,13)*unitscale,'+');
  mms_pa(1).Color = hca.ColorOrder(1,:);
  mms_pa(2).Color = hca.ColorOrder(2,:);
  mms_pa(3).Color = hca.ColorOrder(3,:);
  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLim = [1e-34 2e-27]*unitscale;
  hca.XLim = [1e1 2e4];  
  hca.YLabel.String = 'f_e (s^3 km^{-6})';
  hca.XTick = [1e1 1e2 1e3 1e4];
  %hleg = irf_legend(h,{'0';'90';'180'},[0.98 0.98]);
  hca.Title.String = ['Electron distribution'];
  %hold(h,'on')
  %mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])
  %hold(h,'off')
  irf_legend(hca,{{'-  Fit';'+ Measured';['   ' timeUTC(1:10)];['   ' timeUTC(12:23)]}},[0.05 0.4],'color',[0 0 0])
  
  hold(hca,'on')
  plot(hca,scPot1.resample(time).data*[1 1],hca.YLim)
  hold(hca,'off')
end
if 1 % Measured electron distribution pitchangle circle plot
  hca = h(isub); isub = isub +1;  
  hca.Title.String = ['Electron distribution'];
  unitscale = 1e30; % cm^-6 ->km^-6
  thisPitch = ePitch1(tInd);
  is_photoelectrons = find(EN<scPot1.resample(time).data*1.2);
  [PA,EN] = meshgrid(thisPitch.depend{2},thisPitch.depend{1});
  X = log10(EN).*sind(PA);
  Y = log10(EN).*cosd(PA);
  plotPitch = log10(squeeze(thisPitch.data)*unitscale);
  plotPitch(is_photoelectrons) = NaN;
  pcolor(hca,X,Y,plotPitch)
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hb = colorbar('peer',hca);
  hca.CLim = [-3 0];
  hca.XLabel.String = 'log_{10} E (eV)';
  hca.YLabel.String = 'log_{10} E (eV)';
  hold(hca,'on')
  vdE = units.me*vd(1).^2/2/units.eV;
  hl = plot(hca,log10(vdE),0,'wo');
  hl.MarkerSize = 10;
  hl.LineWidth = 2;
  hold(hca,'off')
end

if 1 % input info
  hca = h(isub); isub = isub + 1;  
  info_str = [ ...
  sprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n)), ...
  sprintf('n = ['), sprintf(' %g',n*1e-6), sprintf('] cc\n'), ...
  sprintf('m = ['), sprintf(' %g',m/me), sprintf('] me\n'), ...
  sprintf('q = ['), sprintf(' %g',q/qe), sprintf('] e\n'), ...
  sprintf('T = ['), sprintf(' %g',T), sprintf('] eV\n'), ...
  sprintf('vt = ['), sprintf(' %.0f',vt*1e-3), sprintf(' ]  km/s\n'), ...
  sprintf('vd = ['), sprintf(' %.0f',vd*1e-3), sprintf(' ]  km/s\n'), ...
  sprintf('wp = ['), sprintf(' %g',wp), sprintf('] Hz\n'), ...
  sprintf('wc = ['), sprintf(' %g',wc), sprintf('] Hz\n'), ...
  sprintf('ro = ['), sprintf(' %g',ro*1e-3), sprintf('] km\n'), ...
  sprintf('Lin = ['), sprintf(' %g',Lin*1e-3), sprintf('] km\n'), ...
  sprintf('Ld = ['), sprintf(' %g',Ld*1e-3), sprintf('] km\n'), ...
  ];
  hca.Visible = 'off';
  text(hca,0,1,info_str,'verticalalignment','top')
end

if 1 % input distributions
  hca = h(isub); isub = isub + 1;
  vmax = max(2*vt + vd);
  vmin = min(-2*vt + vd);
  vvec = linspace(vmin,vmax,1000); 
  
  ftot = 0;
  for isp = 1:nsp
    ftot = ftot + f(vvec,n(isp),vt(isp),vd(isp));
  end
  % plot normalization
  if do_normf
    for isp = 1:nsp
      fnorm(isp) = max(f(vvec,n(isp),vt(isp),vd(isp)));
    end
    fnormtot = max(ftot);
    fmaxtot = 1;
  else
    fnorm = ones(isp,1);
    fnormtot = 1;
    fmaxtot = max(ftot);
  end
  
  % obtained phase velocity
  plot(hca,vphmax*1e-6*[1 1],[0 fmaxtot],'-','linewidth',1.5,'color',0.8+[0 0 0])
  hold(hca,'on')
  

  


  plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  
  colors = [     0    0.4470    0.7410;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0.3010    0.7450    0.9330;...
            0.6350    0.0780    0.1840];
  f_legends = cell(nsp,1);
  for isp = 1:nsp
    %plot(hca,vvec*1e-6,f(vvec,n(isp),vt(isp),vd(isp))/max(f(vvec,n(isp),vt(isp),vd(isp))))    
    %plot(hca,vvec*1e-6,f(vvec,n(isp),vt(isp),vd(isp))/max(f(vvec,n(isp),vt(isp),vd(isp))),'--','linewidth',1.5,'color',colors(isp,:).^0.5)    
    hp = patch(hca,[vvec vvec(end) vvec(1)]*1e-6,[f(vvec,n(isp),vt(isp),vd(isp)) 0 0]/fnorm(isp),colors(isp,:));
    hp.FaceAlpha = 0.3;
    hp.EdgeColor = 'none';
    f_legends{isp} = sprintf('f_%.0f/%g',isp,max(f(vvec,n(isp),vt(isp),vd(isp))));
  end   
  hl = legend(hca,{'v_{ph} @ max w_i',sprintf('f_{tot}/%g',max(ftot)),f_legends{:}},'location','northoutside','box','off');
  hl.Position(3) = 0.75;
  hl.Position(2) = hl.Position(2)+0.08;
  hca.XLabel.String = 'v (10^3 km/s)';
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f'; end
  plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6;  
  hold(hca,'off')
end
if 0 % solution, wr wi, plotyy
  hca = h(isub); isub = isub + 1;    
  if 1
    ax = plotyy(hca,kvec*knorm,wr_store/wnorm,kvec*knorm,wi_store/wnorm);
    %ax(2).YLim(1) = 0;
    ax(2).YGrid = 'on';
    ax(1).XLim = kvec([1 end])*knorm;
    ax(2).XLim = kvec([1 end])*knorm;
    %ax(2).YLim = [-0.03 0.01];
    %ax(2).YTick = [-0.03:0.01:0.01];
    %ax(1).YLim = [0 1.5];
  else
    plot(hca,kvec*knorm,wr_store/wnorm,kvec*knorm,wi_store/wnorm,'linewidth',1.5)
    hca.XLim = [0 max(kvec)*knorm];
  end
  legend(hca,'w_r','w_i', 'location','northeast')
  hca.XLabel.String = sprintf('k%s',knorm_str);
  if wnorm ~= 1
    hca.YLabel.String = sprintf('w/%s',wnorm_str);
  else
    hca.YLabel.String = '2\pi f (2\pi Hz)';
  end
  if wimax>1
    hca.YLim(1) = -wimax/wnorm;
  end
  if 0;%wrmax>1
    hca.YLim(2) = 2*wrmax/wnorm;
  end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
end
if 1 % solution, wr wi, common axis
  hca = h(isub); isub = isub + 1;    
  
  plot(hca,kvec*knorm,wr_store/wnorm,kvec*knorm,wi_store/wnorm,'linewidth',1.5)
  hca.XLim = [0 max(kvec)*knorm];
  hca.YLim = [-0.1 0.3];
  legend(hca,'w_r','w_i', 'location','northeast')
  hca.XLabel.String = sprintf('k%s',knorm_str);
  if wnorm ~= 1
    hca.YLabel.String = sprintf('w/%s',wnorm_str);
  else
    hca.YLabel.String = '2\pi f (2\pi Hz)';
  end
  if 0%wimax>1
    hca.YLim(1) = -wimax/wnorm;
  end
  if 0;%wrmax>1
    hca.YLim(2) = 2*wrmax/wnorm;
  end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  
  ylim = hca.YLim;
  yyaxis(hca,'right')
  hca.YLim = ylim*wp(1)/wp(2);
  hca.YLabel.String = 'w/w_{p2}';
  hca.XLim = [0 0.4];
end
if 1 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec*knorm,vph_store*1e-3,kmax*knorm,vphmax*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end])*knorm;
  legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = sprintf('k%s',knorm_str);
end

%% Electron out flow