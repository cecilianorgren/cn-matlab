units = irf_units;
n0 = 0.035;
T0 = 350; vt0 = sqrt(2*units.e*T0./units.me); % m/s
v0 = 0;
n_psi = 129; min_psi = 0.01; max_psi = 45*T0;
n_vpsi = 100; min_vpsi = 0.01; max_vpsi = 50000;
doLogspace = 1;
if doLogspace
  psi = logspace(log10(min_psi),log10(max_psi),n_psi);
  vpsi = logspace(log10(min_vpsi),log10(max_vpsi),n_vpsi);
else  
  psi = linspace(min_psi,max_psi,n_psi);
  vpsi = linspace(min_vpsi,max_vpsi,n_vpsi);
end

[PSI,VPSI] = ndgrid(psi,vpsi);

n_lb_map = zeros(n_psi,n_vpsi);
n_sep_map = zeros(n_psi,n_vpsi);
v_lb_map = zeros(n_psi,n_vpsi);
v_sep_map = zeros(n_psi,n_vpsi);

for i_psi = 1:n_psi
  fprintf('%g,',i_psi)
  for i_vpsi = 1:n_vpsi
    tmp_psi = psi(i_psi);
    tmp_vpsi = vpsi(i_vpsi);
    [n_lb_tmp,n_sep_tmp,v_lb_tmp,v_sep_tmp] = paper_electron_acceleration.liouville_mapped_nf(n0,T0,v0,tmp_psi,tmp_vpsi);
    n_lb_map(i_psi,i_vpsi) = n_lb_tmp;
    n_sep_map(i_psi,i_vpsi) = n_sep_tmp;
    v_lb_map(i_psi,i_vpsi) = v_lb_tmp;
    v_sep_map(i_psi,i_vpsi) = v_sep_tmp;
  end
end
fprintf('\n')

%% Plot, first
h = setup_subplots(2,1);
isub = 1;

if 0 % colormap
  hca = h(isub); isub = isub + 1;
  pcolor(hca,PSI,VPSI*1e-3,n_sep_map*1e-6)
  shading(hca,'flat');
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi (V)';
  hca.YLabel.String = 'v_{\psi} (10^3 km/s)';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
end
if 1 % contour plot
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI,VPSI*1e-3,n_sep_map*1e-6,0.002:0.002:0.04);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi (V)';
  hca.YLabel.String = 'v_{\psi} (10^3 km/s)';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  plot(2400,17,'r+')
  plot(800,0,'go')
  hold(hca,'off')  
end

if 1 % contour plot, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  plot(2400/T0,17/(vt0*1e-6),'r+')
  plot(800/T0,0,'go')
  hold(hca,'off')  
end

%% Plot, including phi/Tlb statistics
% phi statistics obtained from running paper_electron_acceleration.phi_acc_statistics
% available possibly relevant data
% acc_pot
% tepar_lobe
% te_lobe
% te
% ne_sep
% ne_sep_min

psi_t0 = acc_pot./te_lobe;
nevents = numel(psi_t0);

h = setup_subplots(2,2);
isub = 1;

if 1 % contour plot of n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  %hold(hca,'on')
  %plot(hca,2400/T0,17/(vt0*1e-6),'r+')
  %plot(hca,800/T0,0,'go')
  %hold(hca,'off')  
end
if 1 % contour plot of n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04);
  %clabel(C,hc,'LabelSpacing',100);  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  for ievent = 1:nevents        
    plotx = PSI/T0;            
    ploty = VPSI*1e-3/(vt0*1e-6);    
    plotc = n_sep_map*1e-6;
    
    % use interp1 to get new plotc
    % boundaries for Xq and Yq, new grid points
    minx = psi_t0(ievent);
    maxx = max_psi/T0;
    miny = 0;
    maxy = max_vpsi/vt0;
    
    % interpolate
    %plotc_interp = interp2(plotx,ploty,plotc,Xq,Yq)
    
    [remrow,remcol] = find(plotx<psi_t0(ievent));
    remrow = unique(remrow)';
    remcol = unique(remcol)';
    
    plotx(remrow,:) = [];
    ploty(remrow,:) = [];
    plotc(remrow,:) = [];
    [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep_min(ievent)*[1 1]);
    %clabel(C,hc,'LabelSpacing',100);  
    hc.LineWidth = 1;    
    hc.LineColor = [0 0 0];  
    
%     [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep(ievent)*[1 1]);    
%     hc.LineWidth = 1;
%     hc.LineColor = [0.8 0 0.3];  
  end
  hold(hca,'off')  
  
%   hold(hca,'on')
%   plot(hca,2400/T0,17/(vt0*1e-6),'r+')
%   plot(hca,800/T0,0,'go')
%   hold(hca,'off')  
end
if 1 % contour plot v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),v_sep_map*1e-6,-100:5:100);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{sep,map} (10^3 km/s)';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot v and n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),v_sep_map*1e-6,-100:5:100,'k');
  %clabel(C,hc,'LabelSpacing',100);
  
  hold(hca,'on')
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04,'r');
  %clabel(C,hc,'LabelSpacing',100);
  clabel(C,hc,'Manual');
  hold(hca,'off')
  hold(hca,'on')
  for ievent = 1
    plotx = PSI/T0;            
    ploty = VPSI*1e-3/(vt0*1e-6);    
    plotc = n_sep_map*1e-6;
    
    % use interp1 to get new plotc
    % boundaries for Xq and Yq, new grid points
    minx = psi_t0(ievent);
    maxx = max_psi/T0;
    miny = 0;
    maxy = max_vpsi/vt0;
    
    % interpolate
    %plotc_interp = interp2(plotx,ploty,plotc,Xq,Yq)
    
    [remrow,remcol] = find(plotx<psi_t0(ievent));
    remrow = unique(remrow)';
    remcol = unique(remcol)';
    
    plotx(remrow,:) = [];
    ploty(remrow,:) = [];
    plotc(remrow,:) = [];
    [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep_min(ievent)*[1 1]);
    %clabel(C,hc,'LabelSpacing',100);  
    hc.LineWidth = 1;    
    hc.LineColor = [0 0 0];  
    
%     [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep(ievent)*[1 1]);    
%     hc.LineWidth = 1;
%     hc.LineColor = [0.8 0 0.3];  
  end
  hold(hca,'off')  
  
  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{sep,map} (10^3 km/s)';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
  
  hold(hca,'on')
  plot(hca,2400/T0,17/(vt0*1e-6),'r+')
  plot(hca,800/T0,0,'go')
  hold(hca,'off')  
end
