function [fout,nout] = acc_pot_lio_map(f0_inp,phi_inp)

units = irf_units;
m = units.me;
q = -units.e;
v_dl = -5000*1e3;

% Check and prepare f0
if isa(f0_inp,'TSeries')
  v = f0_inp.depend{1}(1,:)*1e3;
  f_data = mean(f0_inp.data,1);
  
  % Set up f0 by making a gaussian fit to the data
  f0 = fit(v.',f_data.','gauss1');
  
elseif isa(f0_inp,'struct')
  v = f0_inp.v;
  
  
  % Set up f0 based on input parameters
  f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
  f0_str = [f0_str(1:end-1) ';'];
  eval(f0_str)
end

% Check and prepare phi
if isa(phi_inp,'TSeries')
  phi_resamp = phi_inp;%.resample(f0_inp);
  t = phi_resamp.time - phi_resamp.time(1);
  phi = phi_resamp.abs.data;  
end
dv = v(2)-v(1);
[V,PHI] = ndgrid(v,phi);
  
% Initialize variables
F = V*0;
Ffree = V*0;
Ftrap = V*0;
V0 = V*0;
VPH = V*0+v_dl;

% Get free and trapped indices
U = m*(V-VPH).^2/2 + q*PHI;
all_ind = 1:numel(V);
ifree = find(U > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(V-VPH>0);
ibelow = find(V-VPH<0);

V0(iabove) = VPH(iabove) + ((V(iabove)-VPH(iabove)).^2 + 2*q*(PHI(iabove))/m).^0.5; % same as Schamel free streaming
V0(ibelow) = VPH(ibelow) - ((V(ibelow)-VPH(ibelow)).^2 + 2*q*(PHI(ibelow))/m).^0.5;
if v_dl<0
  V0(ibelow) = NaN;
else
  V0(iabove) = NaN;
end
V0(itrap) = NaN;

%U(ibelow) = NaN;

% Get distributions
Ffree(ifree) = f0(V0(ifree));
%Ftrap(itrap) = fsep;
F(ifree) = Ffree(ifree);
%F(itrap) = Ftrap(itrap);

%% Make PDist from obtained F.
fout = f0_inp.clone(phi_inp.time,Ffree');
fout.depend{1} = repmat(fout.depend{1}(1,:),phi_inp.length,1);
nout = irf.ts_scalar(phi_inp.time,nansum(Ffree,1)*dv);
%%

doPlot = 1;
if doPlot
  %%
  h = setup_subplots(3,2);
  
  isub = 1;
  if 1
    hca = h(isub); isub = isub + 1;    
    plot(hca,v*1e-6,f0(v),v*1e-6,f_data)
    hca.XLabel.String = 'v (10^3 km/s)';
    hca.YLabel.String = sprintf('f_e (%s)',f0_inp.units);
    hca.XLim = 30*[-1 1];    
  end
  if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,t,phi)
    hca.XLabel.String = 't-t_0 (s)';
    hca.YLabel.String = 'phi (eV)';
  end
  if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,t,v*1e-6,U/units.eV)
    shading(hca,'flat')
    hca.XLabel.String = 't-t_0 (s)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'U/e (eV)';
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];
  end
  if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,t,v*1e-6,V0)
    shading(hca,'flat')
    hca.XLabel.String = 't-t_0 (s)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'V0 (eV)';
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];
    
  end
  if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,t,v*1e-6,(Ffree))
    shading(hca,'flat')
    hca.XLabel.String = 't-t_0 (s)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'F ()';    
  end
  if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,t,nansum(Ffree,1)*dv*1e-6)
    
    hca.XLabel.String = 't-t_0 (s)';
    hca.YLabel.String = 'n (cm^{-3})';    
    hca.YLim(1) = 0;
  end
end
end