% test flux      
clear h
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels  
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  echannel = 17;
  E_fpi = mean(ePDist1.tlim(tint_phi).depend{1}(:,echannel),1);
  f_485 = squeeze(mean(ePDist1.convertto('s^3/m^6').tlim(tint_phi).data(:,echannel,:,:)));
  ang_polar = ePDist1.tlim(tint_phi).depend{3};
  ang_azim = mean(ePDist1.tlim(tint_phi).depend{2},1);
  d_azim = pi/16;
  d_pol = 1 - cosd(180/16);
  d_vel_vec = linspace(v_edi_minus,v_edi_plus,100);
  d_vel = trapz(d_vel_vec,d_vel_vec.^2); % [v^2dv] = m3/s3
  v_fpi = sqrt(2*units.e*E_fpi./units.me); % m/s
  f_vol = d_vel*d_azim*d_pol;
  flux_fpi_485 = f_485*v_fpi*f_vol; % m-2s-1, (to make into cm-2s-1, multiply with 1e-4)
  flux_scale = 1e6;
  surf(hca,ang_azim,ang_polar,flux_fpi_485'*1e-4/flux_scale); % cm-2s-1
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(flux_scale));
  hca.YLabel.String = 'Polar angle';
  hca.XLabel.String = 'Azimuthal angle';
  hca.Title.String = 'Where electrons are going';
  view(hca,[0 0 1])
  hca.XLim = [0 360];
  hca.YLim = [0 180];
end
if 1
  hca = h(isub); isub = isub + 1;
  [ax,hcb] = mms.plot_skymap(hca,ePDist1,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
end
if 1
  hca = h(isub); isub = isub + 1;
  toplot = ePDist1.d3v('scpot',scPot1.resample(ePDist1));
  [ax,hcb] = mms.plot_skymap(hca,toplot,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
end
if 1
  hca = h(isub); isub = isub + 1;
  toplot = ePDist1.vd3v; toplot.data = toplot.data*1e-6;
  [ax,hcb] = mms.plot_skymap(hca,toplot,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
end

%% Calculate moments from PDist
% particlemoments = mms.psd_moments(pdist,phi,theta,stepTable,energy0,energy1,SCpot,particle,option,option_value)
%particlemoments = mms.psd_moments(ePDist1,ePDist1.depend{2},ePDist1.depend{3},ePDist1.ancillary.esteptable,ePDist1.ancillary.energy0,ePDist1.ancillary.energy1,scPot1,'electron')

pdist = ePDist1;
phi.data = ePDist1.depend{2};
theta.data = ePDist1.depend{3};
stepTable.data = ePDist1.ancillary.esteptable;
energy0 = ePDist1.ancillary.energy0;
energy1 = ePDist1.ancillary.energy1;
SCpot = scPot1;
particle = 'electron';
particlemoments = mms.psd_moments(pdist,phi,theta,stepTable,energy0,energy1,SCpot,particle);

 