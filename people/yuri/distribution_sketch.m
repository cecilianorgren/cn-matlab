% Reduced distribution
units = irf_units;
f_red = @(v,n,T,vd,m) n/((pi)^(1/2)*sqrt(2*units.e*T/m))*exp(-(v-vd).^2./(2*units.e*T/m));

% Wave properties, for marking trapping range
vph = 8000*1e3; % m/s
phi = 50; % eV
vtr = sqrt(2*units.e*phi/units.me); % m/s, for electrons
nv_tr = 100;
vtr_vec = linspace(vph-vtr,vph+vtr,nv_tr);

% Distribution properties
n = [1 1]*1e-6; n_tot = sum(n); % m^-3
T = [200 10]; % eV
vt = sqrt(2*units.e*T/units.me); % m/s, for electrons
vd = [0 14000]*1e3; % m/s
m = [1 1]*units.me; % kg
ns = numel(n); % number of populations

% Constructing f
nv = 500;
v = linspace(-20000,20000,nv)*1e3;
f_sep = zeros(ns,nv); % different populations separate
f_sep_tr = zeros(ns,nv_tr); % different populations separate,
f_sep_vph = zeros(ns,1); 

for is = 1:ns
  f_sep(is,:) = f_red(v,n(is),T(is),vd(is),m(is));
  f_sep_tr(is,:) = f_red(vtr_vec,n(is),T(is),vd(is),m(is));
  f_sep_vph(is,:) = f_red(vph,n(is),T(is),vd(is),m(is));
end
f_tot = sum(f_sep,1); % sum of populations
f_tr_tot = sum(f_sep_tr,1); % sum of populations
f_vph_tot = sum(f_sep_vph,1); % sum of populations


% Plotting
linewidth = 1;
color_tr = [1 0.9 0.4];
facealpha_tr = 0.9;
color_tot = [0 0 0];
color_sep = [0 0 0; 0 0 0];
linestyle_tot = '-';
linestyle_sep = '--';
linestyle_vph = ':';
fontsize_arrow = 10;


hca = subplot(1,1,1);
if 1 % plot f_tot
  h_pl_tot = plot(hca,v*1e-3,f_tot,'linewidth',linewidth,'linestyle',linestyle_tot,'color',color_tot);
end
if 1 % plot f_sep
  hold(hca,'on')
  set(hca,'colororder',color_sep)
  h_pl_sep = plot(hca,v*1e-3,f_sep,'linewidth',linewidth,'linestyle',linestyle_sep);
  hold(hca,'off')
end
if 1 % plot patch corresponding to f_tot to show trapping range
  hold(hca,'on')
  h_tr_patch = patch(hca,[vtr_vec vtr_vec([end 1])]*1e-3,[f_tr_tot 0 0]',color_tr,'linewidth',linewidth,'facealpha',facealpha_tr);    
  hold(hca,'off')
end
if 1 % text and arrows indicating trapping range  
  length = hca.YLim(2)*0.2;
  
  arr_end = [vtr_vec(1)*1e-3 f_tr_tot(1)];  
  arr_start = [arr_end(1) arr_end(2)+length];  
  arrow(arr_start,arr_end)
  text(arr_start(1),arr_start(2),'v_{ph}-v_{tr}  ','horizontalalignment','center','verticalalignment','bottom','fontsize',fontsize_arrow)
  
  arr_end = [vtr_vec(end)*1e-3 f_tr_tot(end)];  
  arr_start = [arr_end(1) arr_end(2)+length];
  arrow(arr_start,arr_end)
  text(arr_start(1),arr_start(2),' v_{ph}+v_{tr}','horizontalalignment','center','verticalalignment','bottom','fontsize',fontsize_arrow)
end
if 1 % text and arrows indicating phase velocity range  
  length = hca.YLim(2)*0.1;
  
  arr_end = [vph(1)*1e-3 0];  
  arr_start = [arr_end(1) arr_end(2) - length];  
  arrow(arr_start,arr_end)
  text(arr_start(1),arr_start(2),'v_{ph}','horizontalalignment','center','verticalalignment','top','fontsize',fontsize_arrow)
  
  hold(hca,'on')
  plot(hca,vph*[1 1]*1e-3,[0 f_vph_tot],'linewidth',linewidth,'linestyle',linestyle_vph)
  hold(hca,'off')
end
if 1 % turn axis off and add zero line at y = 0 instead
  axis(hca,'off')
  % add zero line 
  hold(hca,'on')
  h_v_zero = patch(hca,v*1e-3,v*0,color_tot,'linewidth',linewidth);
  hold(hca,'off')
end
if 1 % add arrow on x-axis, if axis is off
  arr_start = [hca.XLim(end)*0.8 0.0000];  
  arr_end = [hca.XLim(end)*1.00 0.0000];  
  %arrow(arr_start,arr_end,'length',15,'baseangle',90,'tipangle',30)
  arrow(arr_start,arr_end,'baseangle',90,'tipangle',30)
  text(hca.XLim(end)*0.9,0.01,'v_{||}','horizontalalignment','center','verticalalignment','top','fontsize',fontsize_arrow) 
end