% Multi-species streaming instability
% Solver for 1D unmagnetized plasma
% If solution is not as expected, change initial guess: x = 0; 
% i-e instability typically low ~wpi, 
% e-e bump-on-tail instability typically higher ~vd*kvec(ik);
% e-e two-stream instablity typically vbulk*k(1)
% If solution seems ok but out of k range change knorm or k_min, k_max.

%% Set up
units = irf_units;% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps0 = 8.854e-12;
kB = 1.381e-23;
c = 299792458;

vshift = -0.5e7;
% Plasma properties for the different species
B = 10e-9; % not used
n = [0.005 0.008]*1e6;
T = [100 400]; T_K = T*units.eV/kB; % use parallel temperature
m = [1 1]*me;
q = [-1 -1]*qe; 
vd = [-2e7 0] - vshift; % m/s

if 0 % Gradual change
B = 10e-9; % not used
n = [0.005 0.008]*1e6;
T = [400 600]; T_K = T*units.eV/kB; % use parallel temperature
m = [1 1]*me;
q = [-1 -1]*qe; 
vd = [-2.7e7 0.5e7]; % m/s
end


nsp = numel(n); % number of species

% Physical parameters
vt = sqrt(2*qe*T./m); % m/s
wp = sqrt(n.*q.^2./(m*eps0)); % Hz
wc = q*B./m; % not used, although it can affect stability of phase space holes
ro = vt./wc;
Lin = c./wp;
Ld = vt./wp/sqrt(2);

ntot = sum(n); wptot = sqrt(ntot.*q(1).^2./(m(1)*eps0)); % Hz

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
%f = @(v,n,vt,vd) n*(1/pi/vt^2)^(3/2)*exp(-(v-vd).^2/vt.^2);
f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt.^2);
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

%% Dispersion solver, one surface
% Wavenumber vector
nk = 200;
k_min= 0.01; k_max = 2.5;
knorm = min(Ld(1));  % length
knorm_str = sprintf('L_{d%g}',1);
kvec = linspace(k_min,k_max,nk)/knorm;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
x = 000;
for ik = 1:nk  
  xguess = x;
  %xguess = vd(2)*kvec(ik);
  
  af = @(temp) D_streaming(temp,kvec(ik),vt,wp,vd);   
  options = optimset('GradObj','on','display','off','TolFun',1e-16);  
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

%% plot solution, comparing to measured reduced f(vpar)
figure(81)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [10 10 400 1100];

wnorm = sqrt(sum(wp.^2)); wnorm_str = sprintf('w_{p,%s}','tot');

ie = 55;
%time_fe = 
ef = ef1D(ie);
ef.time
do_normf = 0;

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 0 % input info
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

if 1 % input distributions, including observed vph and vtrap
  hca = h(isub); isub = isub + 1;
  vmax = max(2*vt + vd);
  vmin = min(-2*vt + vd);
  vvec = linspace(vmin,vmax,1000); 
    
  colors = mms_colors('matlab');
  
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
  
  % observed distribution
  plot(hca,ef1D.depend{1}(ie,:)*1e-3,ef1D.data(ie,:),'color',colors(3,:)*0,'LineWidth',1.5)
 
  if 1 % add vtrap and vph obs    
    iesw = 4;
    vph = esw_data{9}(iesw);
    c_eval('phi(?) = esw_data{9+?}(iesw);',1:4)
    vtrap = sqrt(2*units.e*phi/units.me)*1e-3; % km/s
    vtrap = max(vtrap);
    time = EpochTT(esw_data{5}{iesw});
    hold(hca,'on')
    plot(hca,-vph*[1 1]*1e-3,hca.YLim,'-','color',1*colors(7,:),'LineWidth',1.5)
    hpatch = patch(hca,[-vph-vtrap -vph-vtrap -vph+vtrap -vph+vtrap]*1e-3,[hca.YLim hca.YLim([2 1])],'c');
    hpatch.FaceAlpha = 0.2;
    hpatch.EdgeColor = hpatch.FaceColor;
    hpatch.EdgeAlpha = hpatch.FaceAlpha;
    hold(hca,'off')    
  end
  hold(hca,'on')
  
  % obtained phase velocity
  plot(hca,(vphmax+vshift)*1e-6*[1 1],hca.YLim,'--','linewidth',1.5,'color',colors(7,:),'LineWidth',1.5)  
  plot(hca,(vvec+vshift)*1e-6,ftot/fnormtot,'-','linewidth',1.5,'color',0*colors(7,:),'linestyle','--')   
  
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
    hp = patch(hca,[vvec vvec(end) vvec(1)]*1e-6+vshift*1e-6,[f(vvec,n(isp),vt(isp),vd(isp)) 0 0]/fnorm(isp),colors(isp,:));
    hp.FaceAlpha = 0.3;
    hp.EdgeColor = 'none';
    f_legends{isp} = sprintf('f_%.0f/%g',isp,max(f(vvec,n(isp),vt(isp),vd(isp))));
  end   
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f (s/m^4)'; end
  %plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6;  
  hca.YLim = [0 1]*1e-3;  
  
  hold(hca,'off')
  %set(gca,'children',flipud(get(gca,'children')))
end
if 0 % input distributions
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
  
  % observed distribution
  colors = mms_colors('matlab');
  plot(hca,ef1D.depend{1}(ie,:)*1e-3-vshift*1e-6,ef1D.data(ie,:),'color',colors(3,:))
  hold(hca,'on')
  
  % obtained phase velocity
  %annotation('arrow',vphmax*1e-6*[1 1],[0 fmaxtot]); 
  %plot(hca,vphmax*1e-6*[1 1],[0 fmaxtot],'-','linewidth',1.5,'color',0.8+[0 0 0])
  plot(hca,vphmax*1e-6*[1 1],hca.YLim,'-','linewidth',1.5,'color',0.8+[0 0 0])
  

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
  %hl = legend(hca,{'v_{ph} @ max w_i',sprintf('f_{tot}/%g',max(ftot)),f_legends{:}},'location','northoutside','box','off');
  %hl.Position(3) = 0.75;
  %hl.Position(2) = hl.Position(2)+0.08;
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f (s/m^4)'; end
  plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6;  
  hca.YLim = [0 1]*1e-3;
  
  hold(hca,'off')
end
if 1 % solution, wr wi, doppler shifted back
  hca = h(isub); isub = isub + 1;    
  if 0
    ax = plot(hca,kvec*knorm,(wr_store+kvec*vshift)/wnorm,kvec*knorm,(wi_store+kvec*vshift)/wnorm);
    %ax(2).YLim(1) = 0;
    ax(2).YGrid = 'on';
    ax(1).XLim = kvec([1 end])*knorm;
    ax(2).XLim = kvec([1 end])*knorm;
    %ax(2).YLim = [-0.03 0.01];
    %ax(2).YTick = [-0.03:0.01:0.01];
    %ax(1).YLim = [0 1.5];
  else
    plot(hca,kvec*knorm,abs((wr_store+kvec*vshift))/wnorm,kvec*knorm,(wi_store+0*kvec*vshift)/wnorm,'linewidth',1.5)
    hca.XLim = [0 max(kvec)*knorm];
    hca.XLim = [0 1];
  end
  legend(hca,'w_r','w_i', 'location','east')
  hca.XLabel.String = sprintf('k%s',knorm_str);
  if wnorm ~= 1
    hca.YLabel.String = sprintf('w/%s',wnorm_str);
  else
    hca.YLabel.String = '2\pi f (2\pi Hz)';
  end
%   if wimax>1
%     hca.YLim(1) = -wimax/wnorm;
%   end
  if 0;%wrmax>1
    hca.YLim(2) = 2*wrmax/wnorm;
  end
  hca.YLim = [hca.YLim(1) 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %text(hc a,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
end
if 1 % solution, wr wi
  hca = h(isub); isub = isub + 1;    
  if 0
    ax = plot(hca,kvec*knorm,wr_store/wnorm,kvec*knorm,wi_store/wnorm);
    %ax(2).YLim(1) = 0;
    ax(2).YGrid = 'on';
    ax(1).XLim = kvec([1 end])*knorm;
    ax(2).XLim = kvec([1 end])*knorm;
    %ax(2).YLim = [-0.03 0.01];
    %ax(2).YTick = [-0.03:0.01:0.01];
    %ax(1).YLim = [0 1.5];
  else
    plot(hca,kvec*knorm,abs(wr_store)/wnorm,kvec*knorm,wi_store/wnorm,'linewidth',1.5)
    hca.XLim = [0 max(kvec)*knorm];
  end
  legend(hca,'w_r','w_i', 'location','east')
  hca.XLabel.String = sprintf('k%s',knorm_str);
  if wnorm ~= 1
    hca.YLabel.String = sprintf('w/%s',wnorm_str);
  else
    hca.YLabel.String = '2\pi f (2\pi Hz)';
  end
%   if wimax>1
%     hca.YLim(1) = -wimax/wnorm;
%   end
  if 0;%wrmax>1
    hca.YLim(2) = 2*wrmax/wnorm;
  end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %text(hc a,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
end
if 0 % solution, wr wi
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
if 0 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec*knorm,vph_store*1e-3,kmax*knorm,vphmax*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end])*knorm;
  legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = sprintf('k%s',knorm_str);
end


%% plot solution
figure(81)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [10 10 400 1100];

wnorm = wp(1); wnorm_str = sprintf('w_{p%g}',1);

do_normf = 0;

nrows = 4;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 0 % input info
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
if 1 % solution, wr wi
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
if 1 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec*knorm,vph_store*1e-3,kmax*knorm,vphmax*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end])*knorm;
  legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = sprintf('k%s',knorm_str);
end
