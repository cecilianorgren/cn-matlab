%% Set up
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
units = irf_units;

% burst interval
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
% shorter time interval for dispersion analysis
Tints = irf.tint('2017-07-06T00:55:39.50Z/2017-07-06T00:55:41.50Z'); % second "good" batch
Tints = irf.tint('2017-07-06T00:54:13.70Z/2017-07-06T00:54:17.00Z'); % first "good" batch

%% Load data
ic = 1:4;
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

% For width amplitude relationship normalization
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('vte? = irf.ts_scalar(ne?.time,sqrt(2*units.eV*(facTe?.xx.data)./units.me)*1e-3); vte?.units = ''km/s'';',ic);
c_eval('vte?par = vte?;',ic);
c_eval('wpe? = irf.ts_scalar(ne?.time,sqrt(units.e^2*(ne?.data*1e6)/units.eps0/units.me)); wpe?.units = ''1/s'';',ic);
c_eval('Lde? = vte?/wpe?/sqrt(2); lDe?.units = ''km'';');
%c_eval('lDe? = irf.ts_scalar(ne?.time,sqrt(units.eps0*units.eV*(facTe?.xx.data)./(ne?.data*1e6)./(ne?.data*1e6)/units.e))*1e3; lDe?.units = ''km'';',ic);

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
gseRav = 0.25*(gseR1 + gseR2 + gseR3 + gseR4);
c_eval('facR? = irf_convert_fac(gseR?-gseRav,gseB?,[0 0 1]);',ic)
c_eval('facB? = irf_convert_fac(gseB?,gseB?,[0 0 1]);',ic)
c_eval('facE? = irf_convert_fac(gseE?,gseB?,[0 0 1]);',ic)
%c_eval('facE?par = irf_convert_fac(gseE?par,gseB?,[0 0 1]);',ic)
 
%% Set up again, after loading data
pathLocalUser = ['/Users/' localuser '/'];
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
dirNameMatlab = sprintf('+mms_%s%s%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
matlabPath = [pathLocalUser '/MATLAB/cn-matlab/' dirNameMatlab '/'];

%% % Load wave phase velocities (the ones obtained semi-manually)
fid = fopen([matlabPath 'esw_properties_redo.txt'],'r');
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

if 1 % load esw_properties_irf_4_v_gui.txt
  fid = fopen([matlabPath 'esw_properties_irf_4_v_gui.txt'],'r');
  data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
  esw_data_ = textscan(fid,[data_format_read]);
  fclose(fid)
  for icell = 1:17%numel(esw_data_)
    esw_data{icell} = cat(1,esw_data{icell},esw_data_{icell});
  end
end

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:17%numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);  
end

tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
c_eval('tsVtrap? = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{9+?}/units.me)*1e-3);',1:4)
c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)

% average
vmax = 0.25*(vmax1 + vmax2 + vmax3 + vmax4); vmax.name = 'av(vph+vtrap)';
vmin = 0.25*(vmin1 + vmin2 + vmin3 + vmin4); vmin.name = 'av(vph-vtrap)';

% max
vmax_data = max([vmax1.data vmax2.data vmax3.data vmax4.data],[],2); 
vmax = irf.ts_scalar(char(esw_data{5}),vmax_data); vmax.name = 'max(vph+vtrap)';
vmin_data = min([vmin1.data vmin2.data vmin3.data vmin4.data],[],2);  
vmin = irf.ts_scalar(char(esw_data{5}),vmin_data); vmin.name = 'min(vph-vtrap)';

%% Find peak to peak scale and transverse instability condition (wb,wce)
t1 = esw_data{3};
t2 = esw_data{4};
lpp = zeros(numel(esw_data{1}),4);
kmax_instability = zeros(numel(esw_data{1}),4);
wb = zeros(numel(esw_data{1}),4);
wg = zeros(numel(esw_data{1}),4);
Te_ = zeros(numel(esw_data{1}),4);
Lde_ = zeros(numel(esw_data{1}),4);
phi_ = [esw_data{10} esw_data{11} esw_data{12} esw_data{13}];

for icell = 1:numel(esw_data{1})
  TT_ = EpochTT([t1{icell}; t2{icell}]);   
  dTT_ = TT_(2) - TT_(1);
  TTc_ = TT_(1) + 0.5*dTT_;
  TT_  = TT_ + 0*dTT_*[+1 -1]; % make time interval slightly shorter, to take away some (one) outliers
  vv_ = esw_data{6}(icell);  
    
  for ic_ = 1:4 
    c_eval('Te_tmp = facTe?.xx.resample(TTc_).data;',ic_);
    c_eval('Lde_tmp = Lde?.resample(TTc_).data;',ic_);
    Te_(icell,ic_) = Te_tmp;
    Lde_(icell,ic_) = Lde_tmp;
    
    c_eval('E_tmp = gseE?par.tlim(TT_);',ic_);
    dt = E_tmp.time(2)-E_tmp.time(1);
    [val_max,ind_max] = max(E_tmp.data);
    [val_min,ind_min] = min(E_tmp.data);
    ii_ = 1:numel(E_tmp.data);
    tt_ = ii_*dt;  
    xx_ = vv_*tt_;   
    lpp(icell,ic_) = abs(abs(ind_max-ind_min)*dt*vv_);
    
    if 0 % lpp(icell,ic_)>100 || lpp(icell,ic_)<10
      %plot(ii_,E_tmp.data,ind_max,val_max,'*',ind_min,val_min,'o')
      plot(xx_,E_tmp.data,xx_(ind_max),val_max,'*',xx_(ind_min),val_min,'o')
      title(sprintf('lpp = %g km',lpp(icell,ic_)))
      pause
    end
    %title(sprintf('ic = %g, icell = %g, lpp = %.0f km',ic_,icell,lpp(icell,ic_)))
    %pause(0.1)    
  end  
  
  units = irf_units;
  B = mean(gseB1.tlim(TT_).abs.data*1e-9);
  kmax_instability(icell,:) = 0.5*pi*units.e*B/units.me.*sqrt(units.me/units.e./abs(phi_(icell,:)));
  wg(icell,:) = units.e*B/units.me;
  wb(icell,:)   = sqrt(4*units.e*abs(phi_(icell,:))/units.me./(0.5*lpp(icell,:)*1e3).^2);
  wb_k(icell,:) = sqrt(4*units.e *abs(phi_(icell,:))/units.me./(lpp(icell,:)*1e3).^2);
end

% remove outlier manually
lpp(17,4) = mean(lpp(17,1:3));
tsLpp = irf.ts_scalar(char(esw_data{5}),lpp);
tsLde = irf.ts_scalar(char(esw_data{5}),Lde_);
tsTe = irf.ts_scalar(char(esw_data{5}),Te_);
% tsPhi = irf.ts_scalar(char(esw_data{5}),phi_); % we already defined above
k = pi./lpp;
kmax_instability = kmax_instability*1e3; % m^-1 -> km^-1
tsK = irf.ts_scalar(char(esw_data{5}),k);
tsK_max_inst = irf.ts_scalar(char(esw_data{5}),kmax_instability);


%%
%scatter(wg(:),wb(:)); axis equal
% Chen 2005, eq. 5 (1D)
% phi and fun_delta are the normalized quantities here
fun_delta = @(phi) sqrt(2.*sqrt(phi)*(4*log(2)-1)./(sqrt(pi).*exp(phi).*(1-erf(sqrt(phi))))); 
  
if 0 % scatterhist, only ESWs within Tints, 1 point for each spacecraft, not normalized  
  % first acceleration channel  
  %h = scatterhist(tsLpp.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:),'Direction','in','Style','stairs');
  h = scatterhist(tsLpp.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:),'Direction','in','Kernel','on','NBins',[10,10]);
  % main plot
  hca = h(1);
  hca.XLabel.String = 'l_{pp} (km)';
  hca.YLabel.String = '\phi (eV)';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:20:200;
  
  drawnow;
  
  % leftmost plot, phi histogram    
  ip = 3;
  h(ip).Position(1) = h(1).Position(1);
  %h(ip).LineWidth = 1.5;
  h(ip).Children(1).LineStyle = 'none';
  h(ip).Children(2).LineStyle = '-.';  
  h(ip).Children(2).LineStyle = '-.';  
  h(ip).Children(2).DisplayName = '\phi';  
  
  % bottom plot, lpp histogram  
  ip = 2;
  h(ip).Position(2) = h(1).Position(2);
  h(ip).Children(1).LineStyle = 'none';
  h(ip).Children(2).DisplayName = 'l_{pp}';
  hleg = legend([h(ip).Children(2) h(3).Children(2)]);
  hleg.Title.String = {'relative','distribution'};
  hleg.Location = 'NorthEast';
  hleg.Position(2) = 0.68;
end
if 0 % scatterhist, only ESWs within Tints, 1 point for each ESW, not normalized  
  % first acceleration channel  
  %h = scatterhist(tsLpp.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:),'Direction','in','Style','stairs');
  h = scatterhist(max(tsLpp.tlim(Tints).data,[],2),max(tsPhi.tlim(Tints).data,[],2),'Direction','in','Kernel','on','NBins',[10,10]);
  %h = scatterhist(mean(tsLpp.tlim(Tints).data,2),mean(tsPhi.tlim(Tints).data,2),'Direction','in','Kernel','on','NBins',[10,10]);
  % main plot
  hca = h(1);
  hca.XLabel.String = 'l_{pp}^{max} (km)';
  hca.YLabel.String = '\phi^{max} (eV)';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:20:200;
  
  drawnow;
  
  % leftmost plot, phi histogram    
  ip = 3;
  h(ip).Position(1) = h(1).Position(1);
  %h(ip).LineWidth = 1.5;
  h(ip).Children(1).LineStyle = 'none';
  h(ip).Children(2).LineStyle = '-.';  
  h(ip).Children(2).LineStyle = '-.';  
  h(ip).Children(2).DisplayName = '\phi';  
  
  % bottom plot, lpp histogram  
  ip = 2;
  h(ip).Position(2) = h(1).Position(2);
  h(ip).Children(1).LineStyle = 'none';
  h(ip).Children(2).DisplayName = 'l_{pp}';
  hleg = legend([h(ip).Children(2) h(3).Children(2)]);
  hleg.Title.String = {'relative','distribution'};
  hleg.Location = 'NorthEast';
  hleg.Position(2) = 0.68;
end
if 1 % scatter plot with errorbar of std for each ESW
  % first acceleration channel: Tints  
  hca = subplot(1,1,1);
  datax = tsLpp.tlim(Tints).data;
  datay = tsPhi.tlim(Tints).data;
  
  ndata = size(datax,1);
  for idata = 1:ndata
    X = mean(datax(idata,:));
    Y = mean(datay(idata,:));
    stdX = std(datax(idata,:));
    stdY = std(datay(idata,:));
    if 0
      YNEG = min(datay(idata,:));
      YPOS = max(datay(idata,:));
      XNEG = min(datax(idata,:));
      XPOS = max(datax(idata,:));
    else
      YNEG = stdY;
      YPOS = stdY;
      XNEG = stdX;
      XPOS = stdX;
    end
    hplot = plot(hca,X,Y,'k.');
    if idata == 1, hold(hca,'on'); end
    errorbar(hca,X,Y,YNEG,YPOS,XNEG,XPOS,'linewidth',1.0,'color',hplot.Color)
  end
  hold(hca,'off')
  
  % main plot
  hca.XLabel.String = 'l_{pp} (km)';
  hca.YLabel.String = '\phi (eV)';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:20:200;
  hca.FontSize = 14;
  
  drawnow;
  hca.Position(1) = 0.2;
  hca.Position(2) = 0.2;
  hca.Position(3) = 0.7;
  hca.Position(4) = 0.7;
end

if 0 % width-amplitude relationship
  
  h = setup_subplots(2,2);
  isub = 1;  
  if 1 % only ESWs within Tints, 1 point for each spacecraft, not normalized
    hca = h(isub); isub = isub + 1;          
    % first acceleration channel
    scatter(hca,tsLpp.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:),'o')
    hold(hca,'on')
    phi_vec = linspace(0,1.4,100);
    plot(hca,1*fun_delta(phi_vec),phi_vec)
    hold(hca,'off')
    hca.XLabel.String = 'l_{pp} (km)';
    hca.YLabel.String = '\phi (eV)';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:20:200;
  end
  if 1 % all ESWs
    hca = h(isub); isub = isub + 1;  
    scatter(hca,lpp(:)/2./Lde_(:),phi_(:)./Te_(:))
    hold(hca,'on')
    % first acceleration channel

    scatter(hca,tsLpp.tlim(Tints).data(:)/2./tsLde.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:)./tsTe.tlim(Tints).data(:),'*')
    phi_vec = linspace(0,1.4,100);
    plot(hca,1*fun_delta(phi_vec),phi_vec)
    hold(hca,'off')
    hca.XLabel.String = '\delta = l_{pp}/2/\lambda_{De}';
    hca.YLabel.String = '\phi = \phi_{max}/T_e';
  end
  if 1 % only ESWs within Tints, 1 point for each spacecraft
    hca = h(isub); isub = isub + 1;          
    % first acceleration channel
    scatter(hca,tsLpp.tlim(Tints).data(:)/2./tsLde.tlim(Tints).data(:),tsPhi.tlim(Tints).data(:)./tsTe.tlim(Tints).data(:),'o')
    hold(hca,'on')
    phi_vec = linspace(0,1.4,100);
    plot(hca,1*fun_delta(phi_vec),phi_vec)
    hold(hca,'off')
    hca.XLabel.String = '\delta = l_{pp}/2/\lambda_{De}';
    hca.YLabel.String = '\phi = \phi_{max}/T_e';
  end
  
end

%% Dispersion analysis
%fhigh = 3000;
%c_eval('facE?_lowpass = facE?.filt(0,fhigh,[],5);')
%c_eval('facE?perp_lowpass = facE?.filt(0,fhigh,[],5);')
load('/Users/cno062/GoogleDrive/Data/polarization_paper_electron_acceleration_perppar.mat')

if 0 
  %%
  Tints_tmp = Tints;
  Tints = Tints +[-2 2];
tic
[xvecs_gse,yvecs_gse,Power_gse] = mms.fk_powerspec4SC('gseE?par.tlim(Tints)','gseR?','gseB?',Tints,'linear',10,'numk',400,'numf',200,'cav',4,'wwidth',2);
[xvecs_perp1,yvecs_perp1,Power_perp1] = mms.fk_powerspec4SC('facE?.tlim(Tints).x','facR?','facB?',Tints,'linear',10,'numk',400,'numf',200,'cav',4,'wwidth',2);
[xvecs_perp2,yvecs_perp2,Power_perp2] = mms.fk_powerspec4SC('facE?.tlim(Tints).y','facR?','facB?',Tints,'linear',10,'numk',400,'numf',200,'cav',4,'wwidth',2);
[xvecs_par,yvecs_par,Power_par] = mms.fk_powerspec4SC('facE?.tlim(Tints).z','facR?','facB?',Tints,'linear',10,'numk',400,'numf',200,'cav',4,'wwidth',2);
toc
Tints = Tints_tmp;
end

%% Compare gse and fac
h(1) = subplot(1,2,1);
h(2) = subplot(1,2,2);

hca = h(1);
pcolor(hca,xvecs_par.kzf,yvecs_par.fkzf,log10(Power_par.Powerkzf)); shading(hca,'flat');

hca = h(2);
pcolor(hca,xvecs_gse.kxf,yvecs_gse.fkxf,log10(Power_gse.Powerkxf)); shading(hca,'flat');

%% Basic plot of dispersion relation
h = setup_subplots(3,3,'vertical');
isub = 1;

knorm = 1e3;

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp1.kxf*knorm,yvecs_perp1.fkxf,log10(Power_perp1.Powerkxf)); shading(hca,'flat');
hca.XLabel.String = 'k_x (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hca.Title.String = 'x = \perp1 ~ z_{GSE}';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp1}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp2.kxf*knorm,yvecs_perp2.fkxf,log10(Power_perp2.Powerkxf)); shading(hca,'flat');
hca.XLabel.String = 'k_x (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp2}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_par.kxf*knorm,yvecs_par.fkxf,log10(Power_par.Powerkxf)); shading(hca,'flat');
hca.XLabel.String = 'k_x (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{||}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp1.kyf*knorm,yvecs_perp1.fkyf,log10(Power_perp1.Powerkyf)); shading(hca,'flat');
hca.XLabel.String = 'k_y (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hca.Title.String = 'y = \perp2 ~ y_{GSE}';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp1}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp2.kyf*knorm,yvecs_perp2.fkyf,log10(Power_perp2.Powerkyf)); shading(hca,'flat');
hca.XLabel.String = 'k_y (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp2}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_par.kyf*knorm,yvecs_par.fkyf,log10(Power_par.Powerkyf)); shading(hca,'flat');
hca.XLabel.String = 'k_y (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{||}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp1.kzf*knorm,yvecs_perp1.fkzf,log10(Power_perp1.Powerkzf)); shading(hca,'flat');
hca.XLabel.String = 'k_z (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hca.Title.String = 'z = || ~ x_{GSE}';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp1}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_perp2.kzf*knorm,yvecs_perp2.fkzf,log10(Power_perp2.Powerkzf)); shading(hca,'flat');
hca.XLabel.String = 'k_z (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{\perp2}';

hca = h(isub); isub = isub + 1;
pcolor(hca,xvecs_par.kzf*knorm,yvecs_par.fkzf,log10(Power_par.Powerkzf)); shading(hca,'flat');
hca.XLabel.String = 'k_z (km^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P_{||}';

hlink = linkprop(h,{'CLim','YLim','XLim'});
hlink.Targets(1).YLim(2) = 1000;
hlink.Targets(1).XLim = [-0.25 0.25];

%% Polarization analysis
tint_polan = Tints + [-1.1 1.1];
ic = 1;
c_eval('Rxyz = gseR?;',ic);
c_eval('Bxyz = gseB?.tlim(tint_polan);',ic);
c_eval('Bscm = gseB?scm.tlim(tint_polan);',ic);
c_eval('Exyz = gseE?.tlim(tint_polan);',ic);


Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
B_SI=Bxyz.abs.data*1e-9;
Wce = e*B_SI/Me;
ecfreq = Wce/2/pi;
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = irf.ts_scalar(Bxyz.time,ecfreq);
ecfreq01 = irf.ts_scalar(Bxyz.time,ecfreq01);
ecfreq05 = irf.ts_scalar(Bxyz.time,ecfreq05);

polarization = irf_ebsp(Exyz,Bscm,Bxyz,Bxyz,Rxyz,[10 4000],'polarization','fac');
%polarization = irf_ebsp(gseE1,gseB1scm,gseB1,gseB1,gseR1,[2 4000],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Bperp = polarization.bb_xxyyzzss(:,:,1)+polarization.bb_xxyyzzss(:,:,2);
Esum = polarization.ee_xxyyzzss(:,:,4);
Eperp = polarization.ee_xxyyzzss(:,:,1)+polarization.ee_xxyyzzss(:,:,2);
Epar = polarization.ee_xxyyzzss(:,:,3);
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);
% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;
vphperp = sqrt(Eperp./Bperp)*1e6;

% Remove points with very low B amplitutes
Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
pfluxz(removepts) = NaN;
vph(removepts) = NaN;
vphperp(removepts) = NaN;

%% Remove electron background
nSecondary = [0.5 1 1.5];
nPhoto = 1;
%[eDist_nobg] = mms.remove_edist_background(ePDist1.);
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''ZeroNaN'',0);',1:numel(nSecondary))


%% Reduced electron distribution
eint = [00 40000];
lowerelim = 000;
vint = [-Inf Inf];
tint_fred = Tints + [-1.1 1.1];
c_eval('eDist = ePDist?.tlim(tint_fred).elim(eint);',ic)
c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
% Remove background
[eDist_nobg] = mms.remove_edist_background(eDist,'tint',tint_fred);
nSecondary = [1.5 2 2.5];
nPhoto = 1;
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''ZeroNaN'',0);',1:numel(nSecondary))


tic; ef1D = eDist_nobg3.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',1000); toc % reduced distribution along B

%% Instability analysis
[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'slow_fit_2'},10,[0.005 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_slow});
[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},344+91i,[0.02 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});
[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},170,[0.013 1 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});

%% Figure: Instability analysis, separate figure;
nrows = 1;
ncols = 2;
h = setup_subplots(nrows,ncols);
h(1).Position(1) = 0.08;
h(1).Position(3) = 0.4;


fontsize = 14;
colors = pic_colors('matlab');
colors =[0 0 0; colors([1:end],:)];
ic = 1;
isub = 1;

% dispersion relation
if 1 % fit for both slow and fast ESWs
  %[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},1000,[0.08 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});
  %[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'slow_fit_2'},10,[0.005 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_slow});  
  
  hca = h(isub); isub = isub + 1;   
  
  % observed slow
  [indmax] = find(wi_1 == max(wi_1));
  kmax = k_1(indmax);
  vphmax = wr_1(indmax)/k_1(indmax);
  tint_disprel_slow_utc = tint_disprel_slow.utc;
  ff = ef1D.tlim(tint_disprel_slow).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_slow).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'*','color',colors(1,:),'displayname',sprintf('electrons:\n%s-%s',tint_disprel_slow_utc(1,12:22),tint_disprel_slow_utc(2,18:22)))
  
  % fit fast
  hold(hca,'on')
  plot(hca,f_1{1}*1e-3,f_1{2},'color',colors(1,:),'linewidth',1.5,'displayname','model input')  
  hold(hca,'off')
  
  % vph for disprel
  hold(hca,'on')
  plot(hca,vphmax*1e-6*[1 1],[0 1.75]*1e-4,'color',colors(1,:),'linewidth',1.5,'displayname','v_{ph} at max growth','linestyle','-.')  
  hold(hca,'off')  
  
  % observed fast
  [indmax] = find(wi_2 == max(wi_2));
  kmax = k_2(indmax);
  vphmax = wr_2(indmax)/k_2(indmax);
  tint_disprel_fast_ = tint_disprel_fast+[0.1 0]+0.05;  
  tint_disprel_fast_utc = tint_disprel_fast_.utc;
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o','color',colors(2,:),'displayname',sprintf('electrons:\n%s-%s',tint_disprel_fast_utc(1,12:22),tint_disprel_fast_utc(2,18:22)))
  
  % fit fast
  plot(hca,f_2{1}*1e-3,f_2{2},'color',colors(2,:),'linewidth',1.5,'displayname','model input')  
  hold(hca,'off')
  
   % vph for disprel
  hold(hca,'on')
  plot(hca,vphmax*1e-6*[1 1],[0 0.61]*1e-4,'color',colors(2,:),'linewidth',1.5,'displayname','v_{ph} at max growth','linestyle','-.')  
  hold(hca,'off')
  
  
  if 1 % ions
    hold(hca,'on')
    ax2 = axes('Position',hca.Position);
    set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
    set(ax2,'YAxisLocation','right');
    set(ax2,'Color','none'); % color of axis
    set(ax2,'XColor','k','YColor','k'); % color of axis lines and numbers
    iscale = 13.5;
    ax2.YLim = [0 hca.YLim(2)]*iscale;
    ax2.YLabel.String = 'f_i (s/m^4)';
    ax2.FontSize = fontsize;
    plot(hca,f_1{1}*1e-3,f_1{3}/iscale,'color',colors(3,:),'linewidth',1.5,'displayname',sprintf('ions:\n%s','model input'))
    hold(hca,'off')
  end
  
  hca.XLim = [-50 50];
  hca.XTick = -100:20:100;
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'f_e (s/m^4)';
  hleg = legend(hca,'location','northwest');
  hleg.Position(1) = 0.09;
  hleg.Position(2) = 0.46;
  hleg.Box = 'off';
  hleg.FontSize = 11;
end
if 0 % fit for fast ESWs
  hca = h(isub); isub = isub + 1;
  
  plot(hca,f_2{1}*1e-3,f_2{2},f_2{1}*1e-3,f_2{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 0 % fit for slow ESWs
  hca = h(isub); isub = isub + 1;  
  irf_pl_mark(h(1),tint_disprel_slow.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_slow,'b'); 
  %[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'slow_fit_2'},10,[0.005 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_slow});  
  plot(hca,f_1{1}*1e-3,f_1{2},f_1{1}*1e-3,f_1{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_slow).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_slow).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 0 % fit for fast ESWs
  hca = h(isub); isub = isub + 1;
  tint_disprel_fast_ = tint_disprel_fast+[0.1 0]+0.05;
  irf_pl_mark(h(1),tint_disprel_fast_.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_fast_,'r')
  %[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},1000,[0.08 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});
  plot(hca,f_2{1}*1e-3,f_2{2},f_2{1}*1e-3,f_2{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 1 % f vs kpar, with all the semi-manual velocities added.
  hca = h(isub); isub = isub + 1;
  kscale = 1e3;
  %pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf)); 
  %pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));  
  hp = pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));
  hp.Annotation.LegendInformation.IconDisplayStyle = 'off';
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'|k| (km^{-1})');
  else
    xlabel(hca,'|k| (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  %colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  k_tint = tsK.tlim(Tints);
  kmax_tint = tsK_max_inst.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs_par.kzf([1 end]);  
  hold(hca,'on')
  % lines with color
  plot(hca,k_1*1e3,wr_1/(2*pi)*1e-3,'color',colors(1,:),'displayname','f_r','linewidth',1.5,'linestyle','-')
  plot(hca,k_1*1e3,20*wi_1/(2*pi)*1e-3,'color',colors(1,:),'displayname','20f_i','linewidth',1.5,'linestyle','--')
  plot(hca,k_2*1e3,wr_2/(2*pi)*1e-3,'color',colors(2,:),'displayname','f_r','linewidth',1.5,'linestyle','-')
  plot(hca,k_2*1e3,1*wi_2/(2*pi)*1e-3,'color',colors(2,:),'displayname','5f_i','linewidth',1.5,'linestyle','--')
  
  % black lines
%   plot(hca,k_1*1e3,wr_1/(2*pi)*1e-3,'k','displayname','f_r')
%   plot(hca,k_1*1e3,20*wi_1/(2*pi)*1e-3,'k--','displayname','10f_i')
%   plot(hca,k_2*1e3,wr_2/(2*pi)*1e-3,'k','displayname','f_r')
%   plot(hca,k_2*1e3,1*wi_2/(2*pi)*1e-3,'k--','displayname','5f_i')
  hold(hca,'off')
  %axis(hca,[0 4e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
  hca.XLim = 0.14*[0 2]; 
  hca.XLim = [0 0.25]; 
  hca.YLim = [0 1];
  colormap(hca,'parula')
  hleg = legend(hca,'Location','East');
  %hleg.Box = 'off';
end
%colormap(h(2),pic_colors('candy'))
%colormap(h(2),'gray')

c_eval('h(?).Position(2) = 0.2; h(?).Position(4) = 0.7;',1:2)
c_eval('ax2.Position(2) = 0.2; ax2.Position(4) = 0.7;',1)
h(2).Position(1) = 0.62;
colormap(h(2),'jet')
c_eval('h(?).FontSize = fontsize;',1:2)

h(2).CLim = [-6 0];
h(2).YLim = [0 0.5];
h(2).XLim = [0 0.18];

%% Plot timeseries and 4SC disprel
npanels = 3;
[h,h2] = initialize_combined_plot(npanels,2,1,0.5,'vertical'); % horizontal

for iip2 = 1:numel(h2)
  h2(iip2).Position(1) = h2(iip2).Position(1) + 0.05;
end

ic = 1;
isub = 0;
zoomy = [];
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % Epar spectra
  isub = isub + 1;
  hca=irf_panel('Eperp spectra');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=0.5*Eperp;
  specrec.f_label='';
  specrec.p_label={'log_{10}E_{\perp}^{2}','mV^2 m^{-2} Hz^{-1}'};
  [~,hcb] = irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
end
if 0 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{tsPhi},'.');
  hca.YLabel.String = {'\phi_{||}','(V)'};  
  hca.YLabel.Interpreter = 'tex';
end

if 0 % v phase
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{tsVphpar},'comp');
  hca.YLabel.String = {'v_{ph}','(km/s)'};  
end
if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');    
  hold(hca,'on')
  
  set(hca,'ColorOrder',mms_colors('111223344'))
  irf_plot(hca,{tsVphpar*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  hca.YLim = sort([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]);
  hold(hca,'off')
  hca.YLabel.String = {'v','(10^3 km/s)'};  
  set(hca,'ColorOrder',mms_colors('122'))
  irf_legend(hca,{'v_{ph}','v_{trap}'},[0.98 0.9],'fontsize',12);    
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('efred');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')

  %set(hca,'ColorOrder',mms_colors('122'))
  %irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  %hca.YLabel.String = 'v_e (km/s)'; 
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
%  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('efred + vtrap');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  vscale = 1e-3;
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))
  irf_plot(hca,{tsVphpar.tlim(Tints)*vscale,vmin.tlim(Tints)*vscale,vmax.tlim(Tints)*vscale},'comp');
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))
  plot(hca,hca.XLim,[0 0],'--','linewidth',1);
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  %hca.YLabel.String = 'v_e (km/s)'; 
  hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(lowerelim,'%g') ' eV'],[0.01 0.99],'color',0*[1 1 1])
  hca.YLim = vscale*sort([max(real(vmax.tlim(Tints).data)) min(real(vmin.tlim(Tints).data))]);
  
  %hca.CLim = [-5 -2.5];
  %hca.CLim = log10(prctile(specrec.p(:),[10 90]));
end

%irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));
irf_zoom(h,'x',Tints + [-1 1]);
irf_zoom(h(zoomy),'y')
%irf_plot_axis_align(h)
%
% dispersion relation
isub = 1;
if 1
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  vph = mms.estimate_phase_speed(Power.Powerkmagf,yvecs.fkmag,xvecs.kmag,100);
  kfit = 2*pi*yvecs.fkmag/vph;
  hold(hca,'on')
  plot(hca,kfit,yvecs.fkmag*1e-3,'k','linewidth',1)
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 1 % f vx |k|, with all the semi-manual velocities added.
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  
  vph_tint = tsVph.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs.kmag([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
  end
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kperp,yvecs.kpar,log10(Power.Powerkperpkpar));
  shading(hca,'flat');
  xlabel(hca,'|k_{\perp}| (m^{-1})');
  ylabel(hca,'k_{||} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_{\perp},k_{||})/P_{max}');
  colormap('jet');
  irf_legend(hca,'(b)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxkxky,yvecs.kykxky,log10(Power.Powerkxky));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kykxky)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kxkxky)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'k_{y} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,k_y)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(c)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxkxkz,yvecs.kzkxkz,log10(Power.Powerkxkz));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kzkxkz)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kxkxkz)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'k_{z} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,k_z)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(d)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kykykz,yvecs.kzkykz,log10(Power.Powerkykz));
  shading(hca,'flat');
  hold(hca,'on')
  plot(hca,[0 0],max(yvecs.kzkykz)*[-1 1],'color','k','Linewidth',1)
  plot(hca,max(xvecs.kykykz)*[-1 1],[0 0],'color','k','Linewidth',1)
  hold(hca,'off')
  xlabel(hca,'k_{y} (m^{-1})');
  ylabel(hca,'k_{z} (m^{-1})');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_y,k_z)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(e)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kxf,yvecs.fkxf*1e-3,log10(Power.Powerkxf));
  shading(hca,'flat');
  xlabel(hca,'k_{x} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_x,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(f)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kyf,yvecs.fkyf*1e-3,log10(Power.Powerkyf));
  shading(hca,'flat');
  xlabel(hca,'k_{y} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_y,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(g)',[0.99 0.99],'color','k','fontsize',14)
end
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kzf,yvecs.fkzf*1e-3,log10(Power.Powerkzf));
  shading(hca,'flat');
  xlabel(hca,'k_{z} (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(k_z,f)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(h)',[0.99 0.99],'color','k','fontsize',14)
end

%% Plot timeseries and 4SC disprel
npanels = 1;
nrows = 1;
ncols = 3;
%h = subplot(nrows,ncols,[1 2]);
%h2(1) = subplot(nrows,ncols,3);

h = subplot('position',[0.07 0.15 0.6 0.8]);
h2 = subplot('position',[0.74 0.15 0.23 0.8]);
%h2(2) = subplot(nrows,ncols,3);

fontsize = 14;
ic = 1;
isub = 1;
zoomy = [];

if 1 % e psd vpar  
  hca = h(isub); isub = isub + 1;
  position = hca.Position;
  vscale = 1e-3;
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;  
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))  
  %pause
  if 0
    irf_plot(hca,{tsVphpar.tlim(Tints)*vscale,vmin.tlim(Tints)*vscale,vmax.tlim(Tints)*vscale},'comp');
  else
    linewidth = 2;
    linemarker = '*';
    linestyle = 'none';
    irf_plot(hca,tsVphpar.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker,'linestyle',linestyle);
    irf_plot(hca,vmin.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
    irf_plot(hca,vmax.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
  end
  
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  if 1 % plot v = 0
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('111'))
    plot(hca,hca.XLim,[0 0],':','linewidth',1);
    hold(hca,'off')
  end
  if 1 % plot vepar
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale,'-');
    hold(hca,'off')
  end
  if 1 % plot vepar +- vth
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+vte1par.resample(gseVe1)*vscale,'--');
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+-1*vte1par.resample(gseVe1)*vscale,'--');
    hold(hca,'off')
  end
  
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';  
  hcb.YLabel.String = 'log_{10} f_e (s/m^4)';
    
  %hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  %hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(lowerelim,'%g') ' eV'],[0.98 0.05],'color',0*[1 1 1],'fontsize',fontsize)
  hca.YLim = vscale*sort([max(real(vmax.tlim(Tints).data)) min(real(vmin.tlim(Tints).data))]);
  hca.YLim = [-40 60];
  
  hca.CLim = [-6 -2.3];
  %hcb.Position(1) = hcb.Position(1) + (position(3) - hca.Position(3));
  %hca.Position = position;
  % h(1).Children
  % to check positions of objects one can do first: legend(hca)
  legend(h(1).Children([3 1 7 5]),{'v_{bulk,||}','v_{bulk,||} \pm v_{te,||}','v_{ph}','v_{ph} \pm v_{tr}'})
end

%irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));
irf_zoom(h(1),'x',Tints + 1*[-1 1]);
irf_zoom(h(zoomy),'y')

% dispersion relation
isub = 1;
if 0
  hca = h2(isub); isub = isub + 1;
  pcolor(hca,xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf)); 
  shading(hca,'flat');
  xlabel(hca,'|k| (m^{-1})');
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');
  irf_legend(hca,'(a)',[0.99 0.99],'color','k','fontsize',14)
  vph = mms.estimate_phase_speed(Power.Powerkmagf,yvecs.fkmag,xvecs.kmag,100);
  kfit = 2*pi*yvecs.fkmag/vph;
  hold(hca,'on')
  plot(hca,kfit,yvecs.fkmag*1e-3,'k','linewidth',1)
  hold(hca,'off')
  axis(hca,[0 3e-4 0 2])
  irf_legend(hca,strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.95 0.95],'color','k','fontsize',14)
  set(gcf,'color','w')
end
if 1 % f vs kpar, with all the semi-manual velocities added.
  hca = h2(isub); isub = isub + 1;
  kscale = 1e3;
  %pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf)); 
  %pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));  
  pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'k_{||} (km^{-1})');
  else
    xlabel(hca,'k_{||} (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  k_tint = tsK.tlim(Tints);
  kmax_tint = tsK_max_inst.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs_par.kzf([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec*kscale,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
    if 1 % k = pi/lpp
      kk_ = min(k_tint.data(ivph,:));
      plot(hca,kk_,kk_/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'ko')
      kk_inst = min(kmax_tint.data(ivph,:));
      %plot(hca,kk_inst,kk_inst/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'k*')
    end
  end
  if 0 % summed up power
    plot(hca,xvecs_par.kzf*kscale,(1/4)*(nansum(Power_par.Powerkzf,1)),'--k','linewidth',0.5)
  end
  hold(hca,'off')
  %axis(hca,[0 4e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
  hca.XLim = 0.14*[-1 1]; 
  hca.YLim = [0 1];
  hleg = legend(hca.Children([end-3 end-4]),{sprintf('%s\n%s\n%s','v_{ph} of','individual','ESWs'),'\lambda\sim\pi/l_{pp}'},'Location','NorthWest');
  hleg = legend(hca.Children([end-3 end-4]),{'v_{ph}','\pi/l_{pp}'},'Location','NorthWest');
  hleg.Box = 'off';
  hleg.Title.String = sprintf('%s\n%s','Parameters of','individual ESWs:');
end
if 0 % f vx |k|, with all the semi-manual velocities added.
  hca = h2(isub); isub = isub + 1;
  kscale = 1e3;
  %pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf)); 
  %pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));
  pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf));
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'|k| (km^{-1})');
  else
    xlabel(hca,'|k| (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  k_tint = tsK.tlim(Tints);
  kmax_tint = tsK_max_inst.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs.kmag([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec*kscale,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
    kk_ = min(k_tint.data(ivph,:));
    plot(hca,kk_,kk_/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'ko')
    kk_inst = min(kmax_tint.data(ivph,:));
    plot(hca,kk_inst,kk_inst/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'k*')
  end
  hold(hca,'off')
  %axis(hca,[0 4e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
  hca.YLim = [0 1000];
end

h.FontSize = fontsize;
h2.FontSize = fontsize;

%% Plot timeseries and 4SC disprel, expanded for instability analysis
npanels = 1;
nrows = 2;
ncols = 3;
%h = subplot(nrows,ncols,[1 2]);
%h2(1) = subplot(nrows,ncols,3);

h(1) = subplot(nrows,ncols,[1 2]);
h(2) = subplot(nrows,ncols,[3]);
h(3) = subplot(nrows,ncols,[4]);
h(4) = subplot(nrows,ncols,[5]);
h(5) = subplot(nrows,ncols,[6]);
%h2(2) = subplot(nrows,ncols,3);

fontsize = 12;
ic = 1;
isub = 1;
zoomy = [];

if 1 % e psd vpar  
  hca = h(isub); isub = isub + 1;
  position = hca.Position;
  vscale = 1e-3;
  specrec = ef1D.specrec('velocity_1D','10^3 km/s');
  specrec.p(specrec.p<1e-7) = NaN;  
  [~,hcb] = irf_spectrogram(hca,specrec);
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('111'))  
  %pause
  if 0
    irf_plot(hca,{tsVphpar.tlim(Tints)*vscale,vmin.tlim(Tints)*vscale,vmax.tlim(Tints)*vscale},'comp');
  else
    linewidth = 2;
    linemarker = '*';
    linestyle = 'none';
    irf_plot(hca,tsVphpar.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker,'linestyle',linestyle);
    irf_plot(hca,vmin.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
    irf_plot(hca,vmax.tlim(Tints)*vscale,'linewidth',linewidth,'marker',linemarker);
  end
  
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  if 1 % plot v = 0
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('111'))
    plot(hca,hca.XLim,[0 0],'--','linewidth',1);
    hold(hca,'off')
  end
  if 1 % plot vepar
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale,'-');
    hold(hca,'off')
  end
  if 1 % plot vepar +- vth
    hold(hca,'on')
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+vte1par.resample(gseVe1)*vscale,'-');
    irf_plot(hca,gseVe1.dot(gseB1.norm.resample(gseVe1))*vscale+-1*vte1par.resample(gseVe1)*vscale,'-');
    hold(hca,'off')
  end
  
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';  
  hcb.YLabel.String = 'log_{10} f_e (s/m^4)';
    
  %hcb.Position(1) = hcb.Position(1) + (h(1).Position(3)-hca.Position(3));
  %hca.Position(3) = h(1).Position(3);
%  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(lowerelim,'%g') ' eV'],[0.98 0.05],'color',0*[1 1 1],'fontsize',fontsize)
  hca.YLim = vscale*sort([max(real(vmax.tlim(Tints).data)) min(real(vmin.tlim(Tints).data))]);
  hca.YLim = [-40 60];
  hca.YLim = [-70 70];
  
  hca.CLim = [-6 -2.5];
  %hcb.Position(1) = hcb.Position(1) + (position(3) - hca.Position(3));
  %hca.Position = position;
end

%irf_zoom(h,'x',irf.tint('2017-07-18T01:40:31.00Z/2017-07-18T01:40:32.50Z'));
irf_zoom(h(1),'x',Tints + 1*[-1 1]);
irf_zoom(h(zoomy),'y')

% dispersion relation

if 1 % f vs kpar, with all the semi-manual velocities added.
  hca = h(isub); isub = isub + 1;
  kscale = 1e3;
  %pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf)); 
  %pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));  
  pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'|k| (km^{-1})');
  else
    xlabel(hca,'|k| (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  k_tint = tsK.tlim(Tints);
  kmax_tint = tsK_max_inst.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs_par.kzf([1 end]);  
  hold(hca,'on')
  for ivph = 1:nvph
    plot(hca,kvec*kscale,kvec/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3,'-k','linewidth',0.5)
    if 1 % k = pi/lpp
      kk_ = min(k_tint.data(ivph,:));
      plot(hca,kk_,kk_/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'ko')
      kk_inst = min(kmax_tint.data(ivph,:));
      plot(hca,kk_inst,kk_inst/(2*pi)*abs(vph_tint.data(ivph)*1e3)*1e-3*1e-3,'k*')
    end
  end
  if 1 % summed up power
    plot(hca,xvecs_par.kzf*kscale,(1/4)*(nansum(Power_par.Powerkzf,1)),'--k','linewidth',0.5)
  end
  hold(hca,'off')
  %axis(hca,[0 4e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
  hca.XLim = 0.14*[-1 1]; 
  hca.YLim = [0 1];
end
colors = pic_colors('matlab');
if 1 % fit for both slow and fast ESWs
  %[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},1000,[0.08 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});
  %[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'slow_fit_2'},10,[0.005 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_slow});  
  
  hca = h(isub); isub = isub + 1;  
  
  tint_disprel_fast_ = tint_disprel_fast+[0.1 0]+0.05;
  irf_pl_mark(h(1),tint_disprel_fast_.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_fast_,'r')  
  irf_pl_mark(h(1),tint_disprel_slow.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_slow,'b');
  
  plot(hca,f_1{1}*1e-3,f_1{2},'color',colors(1,:),'linewidth',1.5)
  if 1 % ions
    hold(hca,'on')
    plot(hca,f_1{1}*1e-3,f_1{3}/15,'color',0*colors(3,:),'linewidth',1.5)
    hold(hca,'off')
  end
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_slow).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_slow).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o','color',colors(1,:))
  hold(hca,'off')
  
  % second fit
  hold(hca,'on')
  plot(hca,f_2{1}*1e-3,f_2{2},'color',colors(2,:),'linewidth',1.5)  
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o','color',colors(2,:))
  hold(hca,'off')
  hca.XLim = [-50 50];
  hca.XTick = -100:20:100;
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'f_{e} (s/m^4)';
  
end
if 0 % fit for fast ESWs
  hca = h(isub); isub = isub + 1;
  
  plot(hca,f_2{1}*1e-3,f_2{2},f_2{1}*1e-3,f_2{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 0 % fit for slow ESWs
  hca = h(isub); isub = isub + 1;  
  irf_pl_mark(h(1),tint_disprel_slow.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_slow,'b'); 
  %[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'slow_fit_2'},10,[0.005 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_slow});  
  plot(hca,f_1{1}*1e-3,f_1{2},f_1{1}*1e-3,f_1{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_slow).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_slow).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 0 % fit for fast ESWs
  hca = h(isub); isub = isub + 1;
  tint_disprel_fast_ = tint_disprel_fast+[0.1 0]+0.05;
  irf_pl_mark(h(1),tint_disprel_fast_.epochUnix,'k')
  irf_pl_mark(h(1),tint_disprel_fast_,'r')
  %[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fast_mod'},1000,[0.08 2 100],{xvecs_par,yvecs_par,Power_par},{ef1D,tint_disprel_fast+[0.1 0]+0.05});
  plot(hca,f_2{1}*1e-3,f_2{2},f_2{1}*1e-3,f_2{3}/10)
  hold(hca,'on')
  ff = ef1D.tlim(tint_disprel_fast_).data;
  ff = mean(ff,1); 
  vv = ef1D.tlim(tint_disprel_fast_).depend{1}(1,:);
  plot(hca,vv*1e-3,ff,'o')
  hold(hca,'off')
end
if 1 % f vs kpar, with all the semi-manual velocities added.
  hca = h(isub); isub = isub + 1;
  kscale = 1e3;
  %pcolor(hca,xvecs_par.kmag*kscale,yvecs_par.fkmag*1e-3,log10(Power_par.Powerkmagf)); 
  %pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));  
  pcolor(hca,xvecs_par.kzf*kscale,yvecs_par.fkzf*1e-3,log10(Power_par.Powerkzf));
  shading(hca,'flat');
  if kscale == 1e3
    xlabel(hca,'|k| (km^{-1})');
  else
    xlabel(hca,'|k| (m^{-1})');
  end
  ylabel(hca,'f (kHz)');
  c=colorbar('peer',hca,'ver');
  ylabel(c,'log_{10} P(f,k)/P_{max}');
  %colormap('jet');  
  
  vph_tint = tsVph.tlim(Tints);
  k_tint = tsK.tlim(Tints);
  kmax_tint = tsK_max_inst.tlim(Tints);
  nvph = vph_tint.length;
  kvec = xvecs_par.kzf([1 end]);  
  hold(hca,'on')
  % lines with color
  plot(hca,k_1*1e3,wr_1/(2*pi)*1e-3,'color',colors(1,:),'displayname','f_r','linewidth',1.5)
  plot(hca,k_1*1e3,20*wi_1/(2*pi)*1e-3,'color',colors(1,:),'displayname','10f_i','linewidth',1.5,'linestyle','--')
  plot(hca,k_2*1e3,wr_2/(2*pi)*1e-3,'color',colors(2,:),'displayname','f_r','linewidth',1.5)
  plot(hca,k_2*1e3,1*wi_2/(2*pi)*1e-3,'color',colors(2,:),'displayname','5f_i','linewidth',1.5,'linestyle','--')
  
  % black lines
%   plot(hca,k_1*1e3,wr_1/(2*pi)*1e-3,'k','displayname','f_r')
%   plot(hca,k_1*1e3,20*wi_1/(2*pi)*1e-3,'k--','displayname','10f_i')
%   plot(hca,k_2*1e3,wr_2/(2*pi)*1e-3,'k','displayname','f_r')
%   plot(hca,k_2*1e3,1*wi_2/(2*pi)*1e-3,'k--','displayname','5f_i')
  hold(hca,'off')
  %axis(hca,[0 4e-4*kscale 0 1])  
  set(gcf,'color','w')
  %axis(hca,'square')
  hca.XLim = 0.14*[0 2]; 
  hca.YLim = [0 1];
  colormap(hca,'parula')
end
colormap(pic_colors('candy'))
irf_zoom(h(1),'x',Tints)
irf_zoom(h(1),'x',tint_fred)
c_eval('h(?).FontSize = fontsize;',1:5)
