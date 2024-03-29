%% Load data to make fit to
localuser = datastore('local','user');
path_dir = ['/Users/' localuser '/GoogleDrive/Research/Michael_Separatrix_EH/f_vs_vpar_72.50/'];
list_of_files = dir([path_dir '*.txt']);
nfiles = numel(list_of_files);
fileids = cell(nfiles,1);
for ifile = 1:nfiles
  fileid = ['f_' list_of_files(ifile).name(7) '_' list_of_files(ifile).name(9:10)];
  fileids{ifile} = fileid;
  data_tmp = load([path_dir list_of_files(ifile).name]);
  eval(sprintf('%s.v = data_tmp(:,1); %s.f = data_tmp(:,2);',fileid,fileid))
end

if 0 % plot distributions
  hca = subplot(1,1,1);
  set(hca,'LineStyleOrder','-|--|:')
  hold(hca,'on')
  for ifile = 1:nfiles
    eval(sprintf('plot(hca,%s.v,%s.f)',fileids{ifile},fileids{ifile}))  
  end
  hca.XLabel.String = 'v_{||}';
  hca.YLabel.String = 'f';
  legend(fileids,'Interpreter','none')
  hold(hca,'off')
end

%% Load all quantities at x = 72.5 as a function of z
path_dir = ['/Users/' localuser '/GoogleDrive/Research/Michael_Separatrix_EH/'];
data_tmp = load([path_dir 'zcut(72.5)_data.txt']);
names_tmp = textread([path_dir 'zcut(72.5)_names.txt'],'%s');
z = data_tmp(:,1);
ne = data_tmp(:,9);
Bx = data_tmp(:,2);
z_picks = [2.43 3.13];
iz = nan(numel(z_picks),1);
for iz_ = 1:numel(z_picks)
  iz(iz_) = find(min(abs(z-z_picks(iz_))) == abs(z-z_picks(iz_)));
end
n_picks = ne(iz); 
B_picks = Bx(iz);
n0 = 0.1;

%figure(83)
%plot(z,ne,z_picks,n_picks,'*')

%% Set up solver, and run through different inputs
units = irf_units;% Physical constants

str_fit_info = '';

dists = [1:2];
ndists = numel(dists);
for idist_ind = 2
  idist = dists(idist_ind);
  [B,n,T,m,q,vd,str_fit_info,x] = get_input(idist);
  nsp = numel(n); % number of species

  % Physical parameters  
  vt = sqrt(2*units.e*T./m); % m/s
  wp = sqrt(n.*q.^2./(m*units.eps0)); % Hz
  wc = q*B./m; % not used, although it can affect stability of phase space holes
  ro = vt./wc;
  Lin = units.c./wp;
  Ld = vt./wp/sqrt(2);

  ntot = sum(n); 
  wptot = sqrt(ntot.*q(1).^2./(m(1)*units.eps0)); % Hz

  disp('----------------------------------------------')
  fprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n))
  fprintf('n = ['); fprintf(' %g',n*1e-6); fprintf('] cc\n');
  fprintf('m = ['); fprintf(' %g',m/units.me); fprintf('] me\n');
  fprintf('q = ['); fprintf(' %g',q/units.e); fprintf('] e\n');
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

  % Dispersion solver, one surface
  % Wavenumber vector
  nk = 100;
  k_min= 0.005; k_max = 0.8;
  knorm = min(Ld(1));  % length
  knorm_str = sprintf('L_{d%g}',1);
  kvec = linspace(k_min,k_max,nk)/knorm;

  wr = nan(1,nk);
  wi = nan(1,nk);
  fval_store = nan(1,nk);

  for ik = 1:nk  
    xguess = x;

    af = @(temp) D_streaming(temp,kvec(ik),vt,wp,vd);   
    options = optimset('GradObj','on','display','off','TolFun',1e-14);  
    [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
    fprintf('x = %g + %gi \n',real(x),imag(x))
    fval_store(ik) = FVAL; 
    wr(ik) = real(x);
    wi(ik) = imag(x);  
  end
  vph = wr./kvec;

  if 0
  rem_ind = nk;350:nk;
  wr_store(rem_ind) = NaN;
  wi_store(rem_ind) = NaN;
  fval_store(rem_ind) = NaN;
  vph(rem_ind) = NaN;
  end
    
  ikmax = find(wi==max(wi),1,'first');
  kmax = kvec(ikmax);
  vphmax = wr(ikmax)/kvec(ikmax);
  wimax = wi(ikmax);
  wrmax = wr(ikmax);

  kmax_all(idist,:) = kmax;
  ikmax_all(idist,:) = ikmax;
  wrmax_all(idist,:) = wrmax;
  wimax_all(idist,:) = wimax;
  vphmax_all(idist,:) = vphmax;
  
  kvec_all(idist,:) = kvec;
  wr_all(idist,:) = wr;
  wi_all(idist,:) = wi;
  vph_all(idist,:) = vph;
end

%% Plot solution, unnormalized
figure(82)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [10 10 1000 1100];

[B,n,T,m,q,vd,str_fit_info,x] = get_input(0);
Va  = B./sqrt(units.mu0*n*units.mp);
%wp = sqrt(n.*q.^2./(m*units.eps0)); % Hz
vt = sqrt(2*units.e*T./m); % m/s

clear h;
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

if 0 % input info
  hca = h(isub); isub = isub + 1;  
  info_str = [ ...
  sprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n)), ...
  sprintf('n = ['), sprintf(' %g',n*1e-6), sprintf('] cc\n'), ...
  sprintf('m = ['), sprintf(' %g',m/units.me), sprintf('] me\n'), ...
  sprintf('q = ['), sprintf(' %g',q/units.e), sprintf('] e\n'), ...
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
if 1 % simulation distributions  
  hca = h(isub); isub = isub + 1;  
  set(hca,'LineStyleOrder','-|--|:')
  hold(hca,'on')
  for ifile = 1:nfiles
    eval(sprintf('plot(hca,%s.v,%s.f)',fileids{ifile},fileids{ifile}))  
  end
  hold(hca,'off')
  hca.XLabel.String = 'v_{||}';
  hca.YLabel.String = 'f';
  legend(hca,fileids,'Interpreter','none','location','northeast')
  hca.XLim = [-20 20];
  hca.Title.String = 'Simulation distributions';
  hca.Box = 'on';
end
if 0 % input distributions,single
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
    %f_legends{isp} = sprintf('f_%.0f/%g',isp,max(f(vvec,n(isp),vt(isp),vd(isp))));
    f_legends{isp} = sprintf('f_%.0f/%g',isp,fnorm(isp));
  end   
  hl = legend(hca,{'v_{ph} @ max w_i',sprintf('f_{tot}/%g',fnormtot),f_legends{:}},'location','northwest','box','off');
  %hl.Position(3) = 0.75;
  %hl.Position(2) = hl.Position(2)+0.08;
  hca.XLabel.String = 'v (10^3 km/s)';
  if do_normf, hca.YLabel.String = 'f/max(f)'; else, hca.YLabel.String = 'f'; end
  %plot(hca,vvec*1e-6,ftot/fnormtot,'-','linewidth',1,'color',[0 0 0])   
  hca.XLim = vvec([1 end])*1e-6; 
  hca.XLim = [-20 20];
  hold(hca,'off')
  hca.Title.String = ['Solver input distribution' str_fit_info];  
end
if 1 % input distributions,all
  hca = h(isub); isub = isub + 1;
  vvec = linspace(-4*vt,4*vt,1000); 
  vvec = linspace(-20*Va,20*Va,1000); 
  ftot(1,:) = get_ftot(1,vvec);
  ftot(2,:) = get_ftot(2,vvec);
  plot(hca,vvec*1e-6,ftot','-','linewidth',1)   
  
  if 1 % phase velocity at max wi  
    hold(hca,'on')
    plotx = vphmax_all*1e-6*[1 1];
    ploty = [1; 1]*[0 hca.YLim(2)];
    colors = get_colors;
    hp = plot(hca,plotx',ploty','-.','linewidth',1.5);
    c_eval('hp(?).Color = colors(?,:).^0.2;',1:numel(vphmax_all))
    hold(hca,'off')   
  end
  colors = get_colors;
  f_legends = cell(nsp,1);  
  %hl = legend(hca,{'v_{ph} @ max w_i'},'location','northwest','box','off');
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = 'f';  
  hca.XLim = vvec([1 end])*1e-6; 
  %hca.XLim = [-20 20];  
  hca.Title.String = 'Solver input distributions';  
end
if 1 % solution, wi
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec_all',wi_all',kmax_all,wimax_all,'*','linewidth',1.5)
  hca.XLabel.String = sprintf('k');
  hca.YLabel.String = '2\pi f_i (2\pi Hz)';
  
  %text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  hca.Title.String = 'Dispersion relation';
end
if 1 % solution, wr
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec_all',wr_all',kmax_all,wrmax_all,'*','linewidth',1.5)
  hca.XLabel.String = sprintf('k');
  hca.YLabel.String = '2\pi f_r (2\pi Hz)';
  
  %text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  hca.Title.String = 'Dispersion relation';
end
if 1 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec_all',vph_all'*1e-3,kmax_all,vphmax_all*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end]);
  %legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = 'k';
end

%% Plot solution, normalized
figure(83)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [10 10 1000 1100];

[B,n,T,m,q,vd,str_fit_info,x] = get_input(0);
Va  = B./sqrt(units.mu0*n*units.mp);
wp = sqrt(n.*q.^2./(m*units.eps0)); % Hz
vt = sqrt(2*units.e*T./m); % m/s
%wc = q*B./m; % not used, although it can affect stability of phase space holes
%ro = vt./wc;
Lin = units.c./wp;
%Ld = vt./wp/sqrt(2);
  
vnorm = Va; vnorm_str = 'v_A';
wnorm = wp; wnorm_str = '\omega_{pe}';
lnorm = Lin; lnorm_str = 'L_{in}'; lnorm_str = 'c/\omega_{pe}'; 

do_normf = 0;

clear h;
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

if 0 % input info
  hca = h(isub); isub = isub + 1;  
  info_str = [ ...
  sprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n)), ...
  sprintf('n = ['), sprintf(' %g',n*1e-6), sprintf('] cc\n'), ...
  sprintf('m = ['), sprintf(' %g',m/units.me), sprintf('] me\n'), ...
  sprintf('q = ['), sprintf(' %g',q/units.e), sprintf('] e\n'), ...
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
if 1 % simulation distributions
  hca = h(isub); isub = isub + 1;  
  set(hca,'LineStyleOrder','-|--|:')
  hold(hca,'on')
  for ifile = 1:nfiles
    eval(sprintf('plot(hca,%s.v,%s.f)',fileids{ifile},fileids{ifile}))  
  end
  hold(hca,'off')
  hca.XLabel.String = 'v_{||}';
  hca.YLabel.String = 'f';
  legend(hca,fileids,'Interpreter','none','location','northeast')
  hca.XLim = [-20 20];
  hca.Title.String = 'Simulation distributions';
  hca.Box = 'on';
end
if 1 % input distributions, all
  hca = h(isub); isub = isub + 1;  
  vvec = linspace(-4*vt,4*vt,1000);  
  ftot(1,:) = get_ftot(1,vvec);
  ftot(2,:) = get_ftot(2,vvec);
  plot(hca,vvec/vnorm,ftot','-','linewidth',1)   
  
  if 1 % phase velocity at max wi  
    hold(hca,'on')
    plotx = vphmax_all*[1 1]/vnorm;
    ploty = [1; 1]*[0 hca.YLim(2)];
    colors = get_colors;
    hp = plot(hca,plotx',ploty','-.','linewidth',1.5);
    c_eval('hp(?).Color = colors(?,:).^0.2;',1:numel(vphmax_all))
    hold(hca,'off')   
  end
  colors = get_colors;
  f_legends = cell(nsp,1);  
  %hl = legend(hca,{'v_{ph} @ max w_i'},'location','northwest','box','off');
  hca.XLabel.String = sprintf('v/%s',vnorm_str);
  hca.YLabel.String = 'f';  
  hca.XLim = vvec([1 end])*1e-6; 
  hca.XLim = [-20 20];  
  hca.Title.String = 'Solver input distributions';  
end
if 1 % solution, wi
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec_all'*lnorm,wi_all'/wnorm,kmax_all*lnorm,wimax_all/wnorm,'*','linewidth',1.5)
  hca.XLabel.String = sprintf('k%s',lnorm_str);  
  hca.YLabel.String = sprintf('w_i/%s',wnorm_str);
  
  %text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  hca.Title.String = 'Dispersion relation';
end
if 1 % solution, wr
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec_all',wr_all',kmax_all,wrmax_all,'*','linewidth',1.5)
  hca.XLabel.String = sprintf('k');
  hca.YLabel.String = '2\pi f_r (2\pi Hz)';
  
  %text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  hca.Title.String = 'Dispersion relation';
end
if 1 % phase velocity
  hca = h(isub); isub = isub + 1;
  ax = plot(hca,kvec_all',vph_all'*1e-3,kmax_all,vphmax_all*1e-3,'*','linewidth',1.5);
  hca.XLim = kvec([1 end]);
  %legend(ax(2),sprintf('v_{ph} @ w_{i,max} = %.0f km/s',vphmax*1e-3),'location','best')
  hca.YLabel.String = 'v_{ph} (km/s)';
  hca.XLabel.String = 'k';
end

%% Plot, 1 solution only
idist_solv = 1; idist_sim = 4; z_str = '2.43';
idist_solv = 2; idist_sim = 7; z_str = '3.13';
figure(84)
fig = gcf;
ud = get(fig);
delete(ud.Children)
fig.Position = [700 200 500 250];

[B,n,T,m,q,vd,str_fit_info,x] = get_input(0);
mp = units.mp; % mass ratio doesn't matter in normalized units since c = 20*vA
Va  = B./sqrt(units.mu0*sum(n)*mp);
c = 20*Va;
wp = sqrt(n.*q.^2./(m*units.eps0)); % Hz
vt = sqrt(2*units.e*T./m); % m/s
Lin = c./wp;
  
vnorm = Va; vnorm_str = 'v_A';
wnorm = wp; wnorm_str = '\omega_{pe}';
lnorm = Lin; lnorm_str = 'L_{in}'; lnorm_str = 'c/\omega_{pe}'; 

do_normf = 0;

clear h;
nrows = 1;
ncols = 2;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

if 0 % input info
  hca = h(isub); isub = isub + 1;  
  info_str = [ ...
  sprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n)), ...
  sprintf('n = ['), sprintf(' %g',n*1e-6), sprintf('] cc\n'), ...
  sprintf('m = ['), sprintf(' %g',m/units.me), sprintf('] me\n'), ...
  sprintf('q = ['), sprintf(' %g',q/units.e), sprintf('] e\n'), ...
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
if 1 % simulation + input distributions
  hca = h(isub); isub = isub + 1;  
  vvec = linspace(-20*Va,20*Va,1000);  
  ftot = get_ftot(idist_solv,vvec);  
  ftot_scale = 0.80e20;
  ftot_scale = 4e21;
  ftot_scale = 1/1e-7;
  eval(sprintf('ftot_scale = 1/((max(ftot)/max(%s.f)));',fileids{idist_sim}))
  set(hca,'LineStyleOrder','-|--|:')    
  eval(sprintf('plot(hca,%s.v,%s.f,vvec/vnorm,ftot*ftot_scale,''LineWidth'',1.5)',fileids{idist_sim},fileids{idist_sim}))    
  hca.XLabel.String = 'v_{||}/v_A';
  hca.YLabel.String = 'f [arb. units]';
  %legend(hca,fileids,'Interpreter','none','location','northeast')
  hca.XLim = [-20 20];
  hca.Title.String = 'Reduced electron distribution';
  %hca.Title.String = {'Reduced','electron distributions'};
  hca.Title.String = {'Reduced electron',sprintf('distribution at z = %s',z_str)};
  hca.Box = 'on';
    
  if 1 % phase velocity at max wi  
    hold(hca,'on')
    plotx = vphmax_all(idist_solv)*[1 1]/vnorm;
    ploty = [0 hca.YLim(2)];
    colors = get_colors;
    hp = plot(hca,plotx',ploty','-','linewidth',1.5,'color',0.8+[0 0 0]);
    %c_eval('hp(?).Color = colors(?,:).^0.2;',1)
    hold(hca,'off')   
  end
  legend(hca,{'Simulation','Fit','v_{ph} at max \omega_i'},'Box','off','location','best')
end
if 1 % solution, wr and wi
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec_all(idist_solv,:)'*lnorm,abs(wr_all(idist_solv,:))'/wnorm,...
           kvec_all(idist_solv,:)'*lnorm,wi_all(idist_solv,:)'/wnorm,...
          'linewidth',1.5)
  hca.XLabel.String = sprintf('k%s',lnorm_str);  
  hca.YLabel.String = ['\omega/' ,wnorm_str];
  hca.YGrid = 'on';
  %irf_legend(hca,{'\omega_r','\omega_i'},[0.01 0.99])
  legend(hca,{'\omega_r','\omega_i'},'Box','off','location','northwest')
  %text(hca,hca.XLim(1),hca.YLim(2),sprintf('w_r @ w_{i,max} = %g\nw_{i,max} = %g',wrmax,wimax),'horizontalalignment','left','verticalalignment','top')
  hca.Title.String = 'Dispersion relation';
end

for ipanel = 1:npanels
  h(ipanel).Position(2) = 0.2;
  h(ipanel).Position(4) = 0.6;
end
h(2).XLim = [0 0.7];
h(2).YLim = [-0.03 0.2];

hl = irf_legend(h(1),'a)',[0.03 0.98],'k'); hl.FontSize = 14;
hl = irf_legend(h(2),'b)',[0.03 0.98],'k'); hl.FontSize = 14;
 
%% Get all normalizations in order
% Normalization
[B0,n0,T0,m0,q0,vd,~,~] = get_input(0);
mp = units.mp;
Va0  = B0./sqrt(units.mu0*n0*mp);
c0 = 20*Va0;
wp0 = sqrt(n0.*q0.^2./(m0*units.eps0)); % Hz
wc0 = units.e*B0/units.me;
vt0 = sqrt(2*units.e*T0./m0); % m/s
Lin0 = c0./wp0;
wp0/wc0;

% Solver input
[B,n,T,m,~,vd,~,~] = get_input(2);
mp = units.mp;
Va  = B./sqrt(units.mu0*n*mp);
Vatot  = B./sqrt(units.mu0*sum(n)*mp);
c = c0;
wp = sqrt(n.*q.^2./(m*units.eps0)); % Hz
wc = units.e*B./m;
wptot = sqrt(sum(n).*q.^2./(m(1)*units.eps0)); % Hz
vt = sqrt(2*units.e*T./m); % m/s
Lin = c./wp;
Lintot = c./wptot;

%%
disp('----------------------------------------------')
fprintf('check quasi-neutrality: sum(qn) = %g\n',sum(q.*n))
fprintf('B0 = %g, B = %g = %g B0\n',B0,B,B/B0)
fprintf('n = ['); fprintf(' %g',n*1e-6); fprintf('] cc = ['); fprintf(' %g',n/n0); fprintf('] n0\n');
fprintf('m = ['); fprintf(' %g',m/units.me); fprintf('] me\n');
fprintf('q = ['); fprintf(' %g',q/units.e); fprintf('] e\n');
fprintf('T = ['); fprintf(' %g',T); fprintf('] eV\n');
fprintf('vt = ['); fprintf(' %.0f',vt*1e-3); fprintf('] km/s = ['); fprintf(' %g',vt/Va0); fprintf('] Va0\n');
fprintf('vd = ['); fprintf(' %.0f',vd*1e-3); fprintf('] km/s = ['); fprintf(' %g',vd/Va0); fprintf('] Va0\n');
fprintf('wp = ['); fprintf(' %g',wp); fprintf('] Hz = ['); fprintf(' %g',wp/wp0); fprintf('] wp0\n');
fprintf('wc = ['); fprintf(' %g',wc); fprintf('] Hz\n');
fprintf('ro = ['); fprintf(' %g',ro*1e-3); fprintf('] km\n');
fprintf('Lin = ['); fprintf(' %g',Lin*1e-3); fprintf('] km = ['); fprintf(' %g',Lin/Lin0); fprintf('] Lin0 (Lin = c/wpe)\n');
fprintf('Ld = ['); fprintf(' %g',Ld*1e-3); fprintf('] km\n');
fprintf('wp0/wc0 = %g, wp/wc = %g\n',wp0/wc0,wp/wc)
fprintf('vA0 = %g m/s, vA = %g m/s\n',Va0,Vatot)
fprintf('vA0 = %g c0, vA = %g c0\n',Va0/c0,Vatot/c0)
disp('----------------------------------------------')
fprintf('wi = ['); fprintf(' %g',wimax_all); fprintf('] Hz = ['); fprintf(' %g',wimax_all/wptot); fprintf('] wp = ['); fprintf(' %g',wimax_all/wp0); fprintf('] wp0\n');
fprintf('wr = ['); fprintf(' %g',wrmax_all); fprintf('] Hz = ['); fprintf(' %g',wrmax_all/wptot); fprintf('] wp = ['); fprintf(' %g',wrmax_all/wp0); fprintf('] wp0\n');
fprintf('k = ['); fprintf(' %g',kmax_all); fprintf('] m-1 = ['); fprintf(' %g',kmax_all*Lintot); fprintf('] Ln-1 = ['); fprintf(' %g',kmax_all*Lin0); fprintf('] Ln0-1 \n');
fprintf('vph = ['); fprintf(' %.0f',vphmax_all*1e-3); fprintf('] km/s = ['); fprintf(' %g',vphmax_all/Va0); fprintf('] Va0\n');
disp('----------------------------------------------')
%% Aid functions
function [B,n,T,m,q,vd,str_fit_info,x] = get_input(index)
  % normalize to simulations and to n0_cc  
  units = irf_units;
  n0 = 0.024e6; 
  B0 = 25e-9; 
  switch index % set up input parameters
    case 0 % normalization
      nref = 1;
      Bref = 1;
      B = B0*Bref; % not used
      n = [0.1]*1e6;
      T = 4*[100]; % use parallel temperature
      m = [1]*units.me;
      q = [-1]*units.e; 
      vd = [0]; % m/s   
      str_fit_info = ''; 
      x = -400; % xguess
      n = nref*n/sum(n);
    case 1 % f_2_43
      nref = 0.1227;
      Bref = 0.6833;
      B = B0*Bref; % not used
      n = [0.08 0.30]*1e6;
      T = 7.5*[70 140]; % use parallel temperature
      m = [1 1]*units.me;
      q = [-1 -1]*units.e; 
      vd = 3.0*[-9000e3 4500e3]; % m/s   
      str_fit_info = ': fit to f at z=2.43'; 
      x = -20; % xguess
      n = nref*n/sum(n);       
    case 2 % f_3_13
      if 0 % these numbers were originally done for fit to a slice of 3D distributiin, NOT 1D reduced, i.e. it was wrong
        nref = 0.0355;
        Bref = 0.8051;
        B = B0*Bref; % not used
        n = [0.01 0.09]*1e6;
        T = 1*8*[15 130]; % use parallel temperature
        m = [1 1]*units.me;
        q = [-1 -1]*units.e; 
        vd = 1*3*[-12000e3 6000e3]; % m/s 
        str_fit_info = ': fit to f at z=3.13';
        x = -10; % xguess
        n = nref*n/sum(n);   
      else
        nref = 0.0355;
        Bref = 0.8051;
        B = B0*Bref; % not used
        n = [0.01 0.01]*1e6;
        T = 1*8*[15 130]; % use parallel temperature
        m = [1 1]*units.me;
        q = [-1 -1]*units.e; 
        vd = 1*3*[-12000e3 6000e3]; % m/s 
        str_fit_info = ': fit to f at z=3.13';
        x = -10; % xguess
        n = nref*n/sum(n);
      end
    case 22 % f_3_13 % copy 2
      nref = 0.0355;
      Bref = 0.6833;
      B = B0*Bref; % not used
      n = [0.01 0.09]*1e6;
      T = 4*[15 130]; % use parallel temperature
      m = [1 1]*units.me;
      q = [-1 -1]*units.e; 
      vd = 2*[-11000e3 5500e3]; % m/s 
      str_fit_info = ': fit to f at z=3.13';
      x = -100; % xguess
      n = nref*n/sum(n); 
    case 3 % f_3_13
      nref = 0.0355;
      Bref = 0.8051;
      B = B0*Bref; % not used
      n = [0.01 0.09]*1e6;
      T = 10*[15 130]; % use parallel temperature
      m = [1 1]*units.me;
      q = [-1 -1]*units.e; 
      vd = 10*[-11000e3 5500e3]; % m/s 
      str_fit_info = ': fit to f at z=3.13';
      x = -100; % xguess
      n = nref*n/sum(n);   
  end
  %[~,n0,~,~,~,~,~,~] = get_input(index);
  %ntot = sum(n);
  %n = n/ntot*n0;
  n = n*n0;
end
function out = get_colors()
  out = [     0    0.4470    0.7410;...
              0.8500    0.3250    0.0980;...
              0.9290    0.6940    0.1250;...
              0.4940    0.1840    0.5560;...
              0.4660    0.6740    0.1880;...
              0.3010    0.7450    0.9330;...
              0.6350    0.0780    0.1840];
  out = [out; out];
end
function out = get_ftot(index,vvec) 
  units = irf_units;
  f = @(v,n,vt,vd) n*(1/pi/vt^2)^(1/2)*exp(-(v-vd).^2/vt.^2);
  [B,n,T,m,q,vd,str_fit_info,x] = get_input(index);
  vt = sqrt(2*units.e*T./m); % m/s
  nsp = numel(n);
  ftot = 0;
  for isp = 1:nsp
    ftot = ftot + f(vvec,n(isp),vt(isp),vd(isp));
  end
  out = ftot;
end
