
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

%% Plots of wave properties
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
if 0 % scatter plot with errorbar of std for each ESW
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
  hca.YLabel.String = '\phi_{max} (V)';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:20:200;
  hca.FontSize = 14;
  
  drawnow;
  hca.Position(1) = 0.2;
  hca.Position(2) = 0.25;
  hca.Position(3) = 0.7;
  hca.Position(4) = 0.7;
end
if 0 % 2 panels: scatter plot with errorbar of std for each ESW, also one panel with normalized values
  % first acceleration channel: Tints  
  h = setup_subplots(1,2);
  hca = h(1);
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
  hca.YLabel.String = '\phi_{max} (V)';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:20:200;
  hca.FontSize = 14;
  
  drawnow;
  %hca.Position(1) = 0.2;
  %hca.Position(2) = 0.25;
  %hca.Position(3) = 0.7;
  %hca.Position(4) = 0.7;
  
  hca = h(2);
  datax_ = tsLpp.tlim(Tints).data;
  datay_ = tsPhi.tlim(Tints).data;
  datax_norm = tsLde.tlim(Tints).data;
  datay_norm = tsTe.tlim(Tints).data;
  
  datax = datax./datax_norm;
  datay = datay./datay_norm;
  
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
  hca.XLabel.String = 'l_{pp}/\lambda_{De} ';
  hca.YLabel.String = 'e\phi_{max}/T_e';
  hca.XLim(1) = 0;
  hca.YLim(1) = 0; 
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = 0:10:200;
  hca.FontSize = 14;
  
  drawnow;
  
  irf_legend(h(1),{'(a)'},[-0.2 1.01],'color',[0 0 0],'fontsize',16)
  irf_legend(h(2),{'(b)'},[-0.2 1.01],'color',[0 0 0],'fontsize',16)
end
if 0 % 3 panels: scatter plot with errorbar of std for each ESW, also one panel with normalized values
  % first acceleration channel: Tints  
  h = setup_subplots(1,3);
  isub = 1;
  if 1 % vph vs lpp
    hca = h(isub); isub = isub + 1;
    datax = tsLpp.tlim(Tints).data;
    datay = tsVphpar.tlim(Tints).data;
    
    disp(sprintf('Lpp [min mean max std] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('Vphpar [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
        
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
      hplot = plot(hca,X,Y*1e-3,'k.');
      if idata == 1, hold(hca,'on'); end
      errorbar(hca,X,Y*1e-3,0*YNEG*1e-3,0*YPOS*1e-3,1*XNEG,1*XPOS,'linewidth',1.0,'color',hplot.Color)
    end
    hold(hca,'off')

    
    % panel formatting
    hca.XLabel.String = 'l_{pp} (km)';
    hca.YLabel.String = 'v_{ph} (10^3 km/s)';
    hca.XLim(1) = 0;
    hca.XLim(2) = 150; 
    hca.YLim(1) = 0; 
    hca.YLim(2) = 40; 
    
    
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:20:200;
    hca.FontSize = 16;
    drawnow;
  end
  if 0 % vtr vs vph
    hca = h(isub); isub = isub + 1;
    datax = tsVphpar.tlim(Tints).data;
    c_eval('datay(:,?) = tsVtrap?.tlim(Tints).data;',1:4)
    datay(2,3) = mean(datay(2,[1 2 4])); % manually override

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
      hplot = plot(hca,X*1e-3,Y*1e-3,'k.');
      if idata == 1, hold(hca,'on'); end
      errorbar(hca,X*1e-3,Y*1e-3,YNEG*1e-3,YPOS*1e-3,0*XNEG*1e-3,0*XPOS*1e-3,'linewidth',1.0,'color',hplot.Color)
    end
    hold(hca,'off')

    hold(hca,'on')
    plot(hca,[0 40],[0 40],'color',[0.7 0.7 0.7])
    hold(hca,'off')

    % panel formatting
    hca.XLabel.String = 'v_{ph} (10^3 km/s)';
    hca.YLabel.String = 'v_{tr} (10^3 km/s)';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.XLim(2) = 40;
    hca.YLim(2) = 40; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:10:200;
    hca.FontSize = 16;
    drawnow;
  end
  if 1 % phi vs lpp
    hca = h(isub); isub = isub + 1;
    datax = tsLpp.tlim(Tints).data;
    datay = tsPhi.tlim(Tints).data;
    
    disp(sprintf('Lpp [min mean max] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('phi [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
    

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
    hca.YLabel.String = '\phi_{max} (V)';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:20:200;
    hca.FontSize = 16;

    drawnow;
  end
  if 1 % ephi/Te vs lpp/lde
    hca = h(isub); isub = isub + 1;
    datax_ = tsLpp.tlim(Tints).data;
    datay_ = tsPhi.tlim(Tints).data;
    datax_norm = tsLde.tlim(Tints).data;
    datay_norm = tsTe.tlim(Tints).data;

    datax = datax./datax_norm;
    datay = datay./datay_norm;
    
    disp(sprintf('Lpp/Lde [min mean max] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('ephi/Te [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
    

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
    hca.XLabel.String = 'l_{pp}/\lambda_{De} ';
    hca.YLabel.String = 'e\phi_{max}/T_e';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:10:200;
    hca.FontSize = 16;
  end
  irf_legend(h(1),{'(a)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
  irf_legend(h(2),{'(b)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
  irf_legend(h(3),{'(c)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
end
if 1 % 3 panels: line-and-dot plot for each ESW, also one panel with normalized values
  % first acceleration channel: Tints  
  h = setup_subplots(1,3);
    set(groot,'defaultAxesColorOrder',pic_colors('matlab'),...
      'defaultAxesLineStyleOrder','+-|o-|v-|s-|p-|h-|d-|^-')  

  isub = 1;
  if 1 % vph vs lpp
    hca = h(isub); isub = isub + 1;
    datax = tsLpp.tlim(Tints).data;
    datay = tsVphpar.tlim(Tints).data;
    
    disp(sprintf('Lpp [min mean max std] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('Vphpar [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
        
    ndata = size(datax,1);
    for idata = 1:ndata
      X = datax(idata,:);
      Y = datay(idata,:)*[1 1 1 1]; % vph 
      hplot = plot(hca,X,Y*1e-3,'k.');
      hplot.MarkerSize = 15;
      if idata == 1, hold(hca,'on'); end
      plot(hca,X,Y*1e-3,'linewidth',1.0,'color',hplot.Color)
    end
    hold(hca,'off')
    
    % panel formatting
    hca.XLabel.String = 'l_{pp} (km)';
    hca.YLabel.String = 'v_{ph} (10^3 km/s)';
    hca.XLim(1) = 0;
    hca.XLim(2) = 150; 
    hca.YLim(1) = 0; 
    hca.YLim(2) = 40;     
    
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:20:200;
    hca.FontSize = 16;
    drawnow;
  end
  if 0 % vtr vs vph
    hca = h(isub); isub = isub + 1;
    datax = tsVphpar.tlim(Tints).data;
    c_eval('datay(:,?) = tsVtrap?.tlim(Tints).data;',1:4)
    datay(2,3) = mean(datay(2,[1 2 4])); % manually override

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
      hplot = plot(hca,X*1e-3,Y*1e-3,'k.');
      if idata == 1, hold(hca,'on'); end
      errorbar(hca,X*1e-3,Y*1e-3,YNEG*1e-3,YPOS*1e-3,0*XNEG*1e-3,0*XPOS*1e-3,'linewidth',1.0,'color',hplot.Color)
    end
    hold(hca,'off')

    hold(hca,'on')
    plot(hca,[0 40],[0 40],'color',[0.7 0.7 0.7])
    hold(hca,'off')

    % panel formatting
    hca.XLabel.String = 'v_{ph} (10^3 km/s)';
    hca.YLabel.String = 'v_{tr} (10^3 km/s)';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.XLim(2) = 40;
    hca.YLim(2) = 40; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:10:200;
    hca.FontSize = 16;
    drawnow;
  end
  if 1 % phi vs lpp
    hca = h(isub); isub = isub + 1;
    datax = tsLpp.tlim(Tints).data;
    datay = tsPhi.tlim(Tints).data;
    
    disp(sprintf('Lpp [min mean max] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('phi [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
    markers = {'','','','','',''};
    ndata = size(datax,1);
    for idata = 1:ndata
      X = (datax(idata,:));
      Y = (datay(idata,:));
      %hplot = plot(hca,X,Y,'.');
      if idata == 1, hold(hca,'on'); end
      plot(hca,X,Y,'linewidth',1.0,'linestyle','none')
    end
    hold(hca,'off')

    % main plot
    hca.XLabel.String = 'l_{pp} (km)';
    hca.YLabel.String = '\phi_{max} (V)';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:20:200;
    hca.FontSize = 16;

    drawnow;
  end
  if 0 % ephi/Te vs lpp/lde
    hca = h(isub); isub = isub + 1;
    datax_ = tsLpp.tlim(Tints).data;
    datay_ = tsPhi.tlim(Tints).data;
    datax_norm = tsLde.tlim(Tints).data;
    datay_norm = tsTe.tlim(Tints).data;

    datax = datax./datax_norm;
    datay = datay./datay_norm;
    
    disp(sprintf('Lpp/Lde [min mean max] = [%g %g %g]',min(mean(datax,2)),mean(mean(datax,2)),max(mean(datax,2))))
    disp(sprintf('ephi/Te [min mean max] = [%g %g %g]',min(mean(datay,2)),mean(mean(datay,2)),max(mean(datay,2))))
    

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
    hca.XLabel.String = 'l_{pp}/\lambda_{De} ';
    hca.YLabel.String = 'e\phi_{max}/T_e';
    hca.XLim(1) = 0;
    hca.YLim(1) = 0; 
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = 0:10:200;
    hca.FontSize = 16;
  end
  irf_legend(h(1),{'(a)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
  irf_legend(h(2),{'(b)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
  irf_legend(h(3),{'(c)'},[-0.2 1.01],'color',[0 0 0],'fontsize',18)
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