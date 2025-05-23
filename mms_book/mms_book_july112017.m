%% Load data
ic = 3;
units = irf_units;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

 %
% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',1:4);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',1:4);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_onecount = iPDist?; iPDist?_onecount.data = (iPDist?_onecount.data./iPDistErr?.data).^2;',ic)

% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',1:4)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gseVexB? = 1e3*cross(1e3*gseVe?,1e-9*gseB?.resample(gseVe?));',1:4) % mV/m

% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

c_eval('gseEVexB? = gseE?.resample(gseVexB?) + gseVexB?''; gseEVexB?.name = ''E+VexB gse'';',1:4)

c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[gseCurlEvexB,gseDivEvexB,gseEvexB,~,~,~] = c_4_j('gseR?brsttime','gseEVexB?');
gseCurlEvexB = gseCurlEvexB*1e-6*1e9; % nT/s

nsm = 5;
[gseCurlEvexB_smo,gseDivEvexB_smo,gseEvexB_smo,~,~,~] = c_4_j(gseR1brsttime,gseR2brsttime,gseR3brsttime,gseR4brsttime,gseEVexB1.smooth(nsm),gseEVexB2.smooth(nsm),gseEVexB3.smooth(nsm),gseEVexB4.smooth(nsm));
gseCurlEvexB_smo = gseCurlEvexB_smo*1e-6*1e9; % nT/s

disp('Done.')
%% Rotate things into new coordinate system 
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9239];
lmn = [L;M;N];


% Rotate data
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)

c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',1:4)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',1:4)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''VExB LMN'';',ic)
c_eval('mvaVexB? = gseVexB?*lmn''; mvaVexB?.name = ''VexB LMN'';',1:4)
c_eval('mvaEVexB? = mvaE?.resample(mvaVexB?)+mvaVexB?''; mvaEVexB?.name = ''E+VexB LMN'';',1:4)
mvaCurlEVexB = gseCurlEvexB*lmn''; gseCurlEvexB.name = 'curl(E+VexB) LMN';
mvaCurlEVexB_smo = gseCurlEvexB_smo*lmn''; gseCurlEvexB_smo.name = 'curl(E+VexB) LMN';
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
%c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
%c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)
%mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)

c_eval('[mvaVe?par,mvaVe?perp] = irf_dec_parperp(mvaB?,mvaVe?);',ic)

c_eval('vAe? = 1e-3*1e-9*gseB?.abs/sqrt(1e6*ne?*units.mu0*units.me)'';',ic)
c_eval('vte? = 1e-3*sqrt(2*units.eV*gseTe?.trace/3/units.me);',ic)




c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)

%% Overview plot, induction

tint_edr = irf.tint('2017-07-11T22:34:00.000Z/2017-07-11T22:34:06.000Z'); %20151112071854
tint_overview = irf.tint('2017-07-11T22:33:29.700Z/2017-07-11T22:34:20.000Z'); %20151112071854
tint_overview = tint_edr + 30*[-1 1]; %20151112071854
tint_dist = irf.tint('2017-07-11T22:34:02.367Z/2017-07-11T22:34:02.397Z');
h = irf_plot(4);

fontsize = 14;
fontsize_leg = 14;
nsm = 1;

if 1 % vi
  hca = irf_panel('vi');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xyz'); 0 0 0; 0.6 0.6 0.6];
  vscale = 1;
  irf_plot(hca,{mvaVi3.x*vscale,mvaVi3.y*vscale,mvaVi3.z*vscale},'comp')
  %hca.YLabel.String = 'v_i (km/s)';
  hca.YLabel.String = {'v_i','(km/s)'};
  irf_legend(hca,{'v_{iL}','v_{iM}','v_{iN}'},[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end
h_del = irf_panel('delete');

if 0 % B
  hca = irf_panel('B');
  hca.LineStyleOrder = {'-','-','-','-','--'};
  hca.ColorOrder = mms_colors('xyz');
  irf_plot(hca,{mvaB3.x,mvaB3.y,mvaB3.z},'comp')
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaB3.z,0}); hp.FaceColor = mms_colors('z'); hp.EdgeColor = mms_colors('z');
   
  hca.NextPlot = "replacechildren";               
  hca.YLabel.String = {'B','(nT)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  hca.ColorOrder = [mms_colors('xyz'); 0.4 0.4 0.4; 0.6 0.6 0.6];
  irf_legend(hca,{'B_{L}','B_{M}','B_{N}'},[0.02 0.98],'fontsize',fontsize_leg)
  %irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end

if 1 % ve
  hca = irf_panel('ve 2');
  hca.LineStyleOrder = {'-','-','-','-','--'};
  hca.ColorOrder = [mms_colors('xyz'); 0.4 0.4 0.4; 0.6 0.6 0.6];
  vscale = 1e-3;
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaVe3.y.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('y'); hp.EdgeColor = mms_colors('y');
  hca.NextPlot = "add";  
  hp = irf_patch(hca,{mvaVe3.x.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('x'); hp.EdgeColor = mms_colors('x');
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaVe3.z.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('z'); hp.EdgeColor = mms_colors('z');
  hca.NextPlot = "add";
  hca.ColorOrder = [0.4 0.4 0.4; 0.6 0.6 0.6];
  irf_plot(hca,{-1*vAe3.smooth(20)*vscale,-1*vte3.smooth(20)*vscale},'comp')  
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  hca.ColorOrder = [mms_colors('xyz'); 0.4 0.4 0.4; 0.6 0.6 0.6];
  %irf_legend(hca,{'v_{eL}','v_{eM}','v_{eN}','v_{Ae}','v_{te}'},[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'v_{eL}','v_{eM}','v_{eN}'},[0.02 0.98],'fontsize',fontsize_leg)
  hca.ColorOrder = [0.4 0.4 0.4; 0.6 0.6 0.6; 0.6 0.6 0.6];
  irf_legend(hca,{'v_{Ae}',' ','v_{te}'}',[0.02 0.53],'fontsize',fontsize_leg)
  %irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end

if 0 % ve
  hca = irf_panel('ve');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xyz'); 0 0 0; 0.6 0.6 0.6];
  vscale = 1e-3;
  irf_plot(hca,{mvaVe3.x*vscale,mvaVe3.y*vscale,mvaVe3.z*vscale,-1*vAe3.smooth(nsm)*vscale,-1*vte3.smooth(nsm)*vscale},'comp')
  %hca.YLabel.String = 'v_e (10^3 km/s)';
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
end
if 1 % E+vexB
  hca = irf_panel('E+vexB');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  hp = irf_patch(hca,{mvaEVexB3.y.smooth(nsm),0},'larger');
  hp.FaceColor = mms_colors('y');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaEVexB3.y.smooth(nsm),0},'smaller');
  hp.FaceColor = mms_colors('x');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'(E+v_e\times B)_M','(mV/m)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  %irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end 
if 0 % curl(E+vexB)
  hca = irf_panel('curl(E+vexB)');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  scale = 1e6;
  hp = irf_patch(hca,{mvaCurlEVexB.smooth(nsm)*scale,0});
  %hp.FaceColor = mms_colors('y');
  %hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  %hp = irf_patch(hca,{mvaCurlEVexB.z.smooth(nsm)*scale,0},'smaller');
  %hp.FaceColor = mms_colors('x');
  %hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'\nabla\times(E+v_e\times B)_N','(10^6 nT/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  irf_legend(hca,{'MMS 1-4'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end
if 0 % curl(E+vexB)
  hca = irf_panel('curl(E+vexB)');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  scale = 1e6;
  hp = irf_patch(hca,{-1*mvaCurlEVexB.z.smooth(nsm)*scale,0},'larger');
  hp.FaceColor = mms_colors('y');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  hp = irf_patch(hca,{-1*mvaCurlEVexB.z.smooth(nsm)*scale,0},'smaller');
  hp.FaceColor = mms_colors('x');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'-\nabla\times(E+v_e\times B)_N','(10^6 nT/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  irf_legend(hca,{'MMS 1-4'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end



c_eval('h(?).FontSize = fontsize;',1:numel(h))
hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
hp = findobj(gcf,'type','patch'); 
hp = hp(end:-1:1); c_eval('hp(?).LineWidth = 1.0;',1:numel(hp))
hp = hp(end:-1:1); c_eval('hp(?).FaceAlpha = 0.3;',1:numel(hp))
irf_zoom(h(1),'x',tint_overview)
irf_zoom(h(3:end),'x',tint_edr)
irf_zoom(h,'y')
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))
c_eval('h(?).LineWidth = 1.0;',1:numel(h))
delete(h_del)
%c_eval('h(?).Position(2) = h(?).Position(2)+0.12;',3:numel(h))
c_eval('h(?).Position(2) = h(?).Position(2)+0.15;',3:numel(h))
[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(1),h(3));
h(1).XLabel = [];

irf_pl_mark(h(1),tint_edr(1),'k')
irf_pl_mark(h(1),tint_edr(2),'k')

irf_pl_mark(h(3:end),tint_dist','k')

c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',[1 3:numel(h)])
irf_plot_axis_align(h([1 3:end]))

%% Overview plot, induction + rotation of ve

tint_edr = irf.tint('2017-07-11T22:34:00.000Z/2017-07-11T22:34:06.000Z'); %20151112071854
tint_overview = irf.tint('2017-07-11T22:33:29.700Z/2017-07-11T22:34:20.000Z'); %20151112071854
tint_overview = tint_edr + 30*[-1 1]; %20151112071854
%h = irf_plot(5);
[h,h2] = initialize_combined_plot('leftright',5,2,1,0.6,'horizontal');

fontsize = 14;
fontsize_leg = 12;
nsm = 5;

if 1 % vi
  hca = irf_panel('vi');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xyz'); 0 0 0; 0.6 0.6 0.6];
  vscale = 1;
  irf_plot(hca,{mvaVi3.x*vscale,mvaVi3.y*vscale,mvaVi3.z*vscale},'comp')
  %hca.YLabel.String = 'v_i (km/s)';
  hca.YLabel.String = {'v_i','(km/s)'};
  irf_legend(hca,{'v_{iL}','v_{iM}','v_{iN}'},[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end
h_del = irf_panel('delete');

if 1 % ve
  hca = irf_panel('ve 2');
  hca.LineStyleOrder = {'-','-','-','-','--'};
  hca.ColorOrder = [mms_colors('xyz'); 0.4 0.4 0.4; 0.6 0.6 0.6];
  vscale = 1e-3;
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaVe3.y.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('y'); hp.EdgeColor = mms_colors('y');
  hca.NextPlot = "add";  
  hp = irf_patch(hca,{mvaVe3.x.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('x'); hp.EdgeColor = mms_colors('x');
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaVe3.z.smooth(nsm)*vscale,0}); hp.FaceColor = mms_colors('z'); hp.EdgeColor = mms_colors('z');
  hca.NextPlot = "add";
  hca.ColorOrder = [0.4 0.4 0.4; 0.6 0.6 0.6];
  irf_plot(hca,{-1*vAe3.smooth(20)*vscale,-1*vte3.smooth(20)*vscale},'comp')  
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  hca.ColorOrder = [mms_colors('xyz'); 0.4 0.4 0.4; 0.6 0.6 0.6];
  irf_legend(hca,{'v_{eL}','v_{eM}','v_{eN}','v_{Ae}','v_{te}'},[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end
if 0 % ve
  hca = irf_panel('ve');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xyz'); 0 0 0; 0.6 0.6 0.6];
  vscale = 1e-3;
  irf_plot(hca,{mvaVe3.x*vscale,mvaVe3.y*vscale,mvaVe3.z*vscale,-1*vAe3.smooth(nsm)*vscale,-1*vte3.smooth(nsm)*vscale},'comp')
  %hca.YLabel.String = 'v_e (10^3 km/s)';
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
end
if 1 % E+vexB
  hca = irf_panel('E+vexB');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  hp = irf_patch(hca,{mvaEVexB3.y.smooth(nsm),0},'larger');
  hp.FaceColor = mms_colors('y');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  hp = irf_patch(hca,{mvaEVexB3.y.smooth(nsm),0},'smaller');
  hp.FaceColor = mms_colors('x');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'(E+v_e\times B)_M','(mV/m)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  irf_legend(hca,{'MMS 3'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end 
if 0 % curl(E+vexB)
  hca = irf_panel('curl(E+vexB)');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  scale = 1e6;
  hp = irf_patch(hca,{mvaCurlEVexB.smooth(nsm)*scale,0});
  %hp.FaceColor = mms_colors('y');
  %hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  %hp = irf_patch(hca,{mvaCurlEVexB.z.smooth(nsm)*scale,0},'smaller');
  %hp.FaceColor = mms_colors('x');
  %hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'\nabla\times(E+v_e\times B)_N','(10^6 nT/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  irf_legend(hca,{'MMS 1-4'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end
if 1 % curl(E+vexB)
  hca = irf_panel('curl(E+vexB)');
  hca.LineStyleOrder = {'-','-','-','-','-'};
  hca.ColorOrder = [mms_colors('xxzz'); 0 0 0; 0.6 0.6 0.6];
  %hca.ColorOrder = [mms_colors('1'); 0.4 0.4 0.4; 0.6 0.6 0.6];  
  scale = 1e6;
  %hp = irf_patch(hca,{-1*mvaCurlEVexB_smo.z.smooth(nsm)*scale,0},'larger');
  hp = irf_patch(hca,{-1*mvaCurlEVexB.z.smooth(nsm)*scale,0},'larger');
  hp.FaceColor = mms_colors('y');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.NextPlot = "add";
  %hp = irf_patch(hca,{-1*mvaCurlEVexB_smo.z.smooth(nsm)*scale,0},'smaller');
  hp = irf_patch(hca,{-1*mvaCurlEVexB.z.smooth(nsm)*scale,0},'smaller');
  hp.FaceColor = mms_colors('x');
  hp.EdgeColor = [0.6 0.6 0.6];
  hca.YLabel.String = {'-\nabla\times(E+v_e\times B)_N','(10^6 nT/s)'};
  hca.YLabel.Interpreter = 'tex';
  hca.NextPlot = "replacechildren";
  irf_legend(hca,{'MMS 1-4'},[0.98 0.98],'fontsize',fontsize_leg,'color','k')  
end



c_eval('h(?).FontSize = fontsize;',1:numel(h))
hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
hp = findobj(gcf,'type','patch'); 
hp = hp(end:-1:1); c_eval('hp(?).LineWidth = 1.0;',1:numel(hp))
hp = hp(end:-1:1); c_eval('hp(?).FaceAlpha = 0.3;',1:numel(hp))
irf_zoom(h(1),'x',tint_overview)
irf_zoom(h(3:end),'x',tint_edr)
irf_zoom(h,'y')
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))
c_eval('h(?).LineWidth = 1.0;',1:numel(h))
c_eval('h(?).Position(1) = 0.1;',1:numel(h))
delete(h_del)
c_eval('h(?).Position(2) = h(?).Position(2)+0.12;',3:numel(h))
[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(1),h(3));
h(1).XLabel = [];

irf_pl_mark(h(1),tint_edr(1),'k')
irf_pl_mark(h(1),tint_edr(2),'k')

c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',[1 3:numel(h)])

%%
isub = 1; 
hca = h2(isub); isub = isub + 1;

tint_hodo = tint_edr + 3*[-1 1];
for ic = 1:4
  c_eval('ve = mvaVe?.tlim(tint_hodo).smooth(5)*1e-3;',ic)
  plot(hca,ve.x.data,ve.y.data,'color',mms_colors(num2str(ic)))
  if ic == 1; hca.NextPlot = 'add'; end
end
hca.NextPlot = 'replacechildren';
hca.XLim = 1e-3*[-10000 10000]*0.99;
hca.YLim = 1e-3*[-18000 2000]*0.99;
hca.XTick = -20:2:20;
axis(hca,'square');
hca.Box = 'on';
hca.FontSize = fontsize;
hca.XLabel.String = 'v_{eL} (10^3 km/s)';
hca.YLabel.String = 'v_{eM} (10^3 km/s)';
hca.ColorOrder = mms_colors('1234');
irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'},[0.02 1.01],'fontsize',fontsize_leg);

hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))

hca.NextPlot = 'add';
rx = 0;
ry = -9;
rr = 7.5;
%plot(hca,rx+rr*cosd(0:360),ry+rr*sind(0:360),'--','color',[0.5 0.5 0.5])
angles = 135:405;
vecx = [rx+rr*cosd(angles) rx];
vecy = [ry+rr*sind(angles) ry];
%patch(hca,vecx,vecy,[0.5 0.5 0.5],'facealpha',0.05)
plot(hca,hca.XLim,[0 0],'color',[0.5 0.5 0.5])
plot(hca,[0 0],hca.YLim,'color',[0.5 0.5 0.5])
hca.NextPlot = 'replacechildren';

hca.NextPlot = 'add';
rx = 0;
ry = 0;
rr = 17;
%plot(hca,rx+rr*cosd(0:360),ry+rr*sind(0:360),'--','color',[0.5 0.5 0.5])
angles = 135:405;
vecx = [rx+rr*cosd(angles)];
vecy = [ry+rr*sind(angles)];
%patch(hca,vecx,vecy,[0.5 0.5 0.5],'facealpha',0.05)
plot(hca,vecx,vecy,'color',[0.5 0.5 0.5])
hca.NextPlot = 'replacechildren';


h2 = h2(1);
c_eval('h2(?).FontSize = fontsize;',1:numel(h2))


%% Make reduced distributions
disp('Preparing reduced distributions.')
vint = [-Inf Inf];

if 0 % ion LMN, entire
  c_eval('if1DL? = iPDist?.reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DM? = iPDist?.reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DN? = iPDist?.reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end
if 0 % electron LMN, tlim
  tint_efred = irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:40.00Z');
  c_eval('ef1DL? = ePDist?.tlim(tint_efred).reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DM? = ePDist?.tlim(tint_efred).reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DN? = ePDist?.tlim(tint_efred).reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end
if 1 % ion LMN, ielim
  tint_ifred =  irf.tint('2017-07-11T22:32:00.00Z/2017-07-11T22:35:20.00Z');
  ielim = [1000 40000];
  c_eval('if1DL?_elim = iPDist?.elim(ielim).tlim(tint_ifred).reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DM?_elim = iPDist?.elim(ielim).tlim(tint_ifred).reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DN?_elim = iPDist?.elim(ielim).tlim(tint_ifred).reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end
if 0 % ion LMN, nobg   
  tint_ifred =  irf.tint('2017-07-11T22:32:00.00Z/2017-07-11T22:35:20.00Z');  
  c_eval('if1DL?_nobg = iPDist?_nobg.tlim(tint_ifred).reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DM?_nobg = iPDist?_nobg.tlim(tint_ifred).reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('if1DN?_nobg = iPDist?_nobg.tlim(tint_ifred).reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end
if 0 % electron LMN, elim
  tint_efred =  irf.tint('2017-07-11T22:33:40.00Z/2017-07-11T22:34:30.00Z');
  elim = [100 40000];
  c_eval('ef1DL?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DM?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DN?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end

%% Estimating the ion components
tint_fit =  irf.tint('2017-07-11T22:33:20.00Z/2017-07-11T22:34:30.00Z');
tint_fit =  irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:20.00Z');
vdf = if1DN1_elim.tlim(tint_fit);
%vdf = if1DN1_nobg.tlim(tint_fit);
vdf = vdf(1:1:vdf.length);
%vdf = vdf(fix(vdf.length/2));
%vdf = vdf([1 vdf.length]);
nPop = 3; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  0,      1000e3];

% nPop = 2; % Two populations looks ok, but beams are warmer and slower to cover region around v=0
% X0 = [0.01e6, -1000e3, 500e3,...
%       0.01e6,  1000e3, 500e3];

if 1
nPop = 4; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  -1000e3,      1000e3,...
      0.01e6,  1000e3,      1000e3];
end
tic; [fitdata_,ts_] = funFitVDF(vdf(100),'nPop',nPop,'plot',1,'guessprevious',0,'X0',X0,'weight',repmat([0 1 1],1,nPop)); toc;

%% Plot fit
h = irf_plot(7);
hca = irf_panel('vdf_obs');
irf_spectrogram(hca,vdf.specrec('velocity'),'lin')
hca = irf_panel('vdf_fit');
irf_spectrogram(hca,ts.f.specrec('velocity'),'lin')
hold(hca,'on'); irf_plot(hca,ts.vd,':'); hold(hca,'off'); 
hca = irf_panel('n_fit');
irf_plot(hca,ts.n)
hca = irf_panel('vd_fit');
irf_plot(hca,ts.vd)
hca = irf_panel('T_fit');
irf_plot(hca,ts.T)
hca = irf_panel('T*n_fit');
irf_plot(hca,ts.T.*ts.n)
hca = irf_panel('cost function');
irf_plot(hca,ts.cf)

irf_plot_axis_align(h)
irf_zoom(h,'x',vdf.time([1 vdf.length]))
hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
colormap(pic_colors('candy4'))

%% 2D fit
%tint_vdf =  irf_time('2017-07-11T22:33:50.00Z','utc>EpochTT');
%tint_fit =  irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:00.00Z');
tint_fit =  irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:00.00Z');
tint_vdf =  irf_time('2017-07-11T22:34:00.00Z','utc>EpochTT');
tint_vdf = tint_fit([1 2]);

ielim = [1000 40000];
vg = -2500:100:2500; % higher grid resolution slows down fit
c_eval('vdf = iPDist?.elim(ielim).tlim(tint_vdf).reduce(''2D'',M,N,''vint'',vint,''scpot'',scPot?,''vg'',vg);',ic)
vdf = vdf([1]);
tint_fit =  irf.tint('2017-07-11T22:33:20.00Z/2017-07-11T22:34:30.00Z');

nPop = 2;
X0 = [0.01e6, -1000e3, -500e3, 500e3, 500e3,...
      0.01e6,  1000e3, -500e3, 500e3, 500e3];
nPop = 3;
X0 = [0.01e6, 500e3, -200e3, 200e3, 200e3,... % initial xy/MN switched..?
      0.01e6,  -500e3, -200e3, 200e3, 200e3,...
      0.02e6,  000e3, 000e3, 2000e3, 2000e3];

nPop = 4;
X0 = [0.01e6, 500e3, -200e3, 200e3, 200e3,... % initial xy/MN switched..?
      0.01e6,  -500e3, -200e3, 200e3, 200e3,...
      0.02e6,  500e3, 000e3, 2000e3, 2000e3,...
      0.02e6,  -500e3, 000e3, 2000e3, 2000e3];    

% nPop = 4;
% X0 = [0.01e6, -500e3, -500e3, 800e3, 800e3,...
%       0.01e6,  -500e3, 500e3, 800e3, 800e3,...
%       0.01e6,  000e3, -500e3, 1000e3, 1000e3,...
%       0.01e6,  000e3, 500e3, 1000e3, 1000e3];
[fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',0,'guessprevious',0,'X0',X0,'maxFunEvals',2000);
%%
iplot = 1;
h = setup_subplots(1,2+nPop); isub = 1;

hca = h(isub); isub = isub + 1;
imagesc(hca,vdf.depend{1}(iplot,:),vdf.depend{1}(iplot,:),squeeze(vdf.data(iplot,:,:)))

hca = h(isub); isub = isub + 1;
imagesc(hca,ts.f.depend{1}(iplot,:),ts.f.depend{1}(iplot,:),squeeze(ts.f.data(iplot,:,:)))

for iPop = 1:nPop  
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,ts.f_sep{iPop}.depend{1}(iplot,:),ts.f_sep{iPop}.depend{1}(iplot,:),squeeze(ts.f_sep{iPop}.data(iplot,:,:)))
end

hlinks = linkprop(h,{'CLim'});
colormap(pic_colors('thermal'))
c_eval('axis(h(?),''square'');',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).XLabel.String = ''v_N''; h(?).YLabel.String = ''v_M'';',1:numel(h))

%% Check if 1D and 2D reduced makes sense, looks good
tint_vdf =  irf_time('2017-07-11T22:34:00.00Z','utc>EpochTT');
ielim = [1000 40000];
vg = -2500:100:2500;

% vdf_2D_LM = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('2D',L,M,'vint',vint,'vg',vg);
% vdf_1D_L = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('1D',L,'vint',vint,'vg',vg);
% vdf_1D_M = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('1D',M,'vint',vint,'vg',vg);
vdf_2D_MN = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('2D',M,N,'vint',vint,'vg',vg);
vdf_1D_M = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('1D',M,'vint',vint,'vg',vg);
vdf_1D_N = iPDist1.elim(ielim).tlim(tint_vdf+0.076*[-1 1]).reduce('1D',N,'vint',vint,'vg',vg);

ielim = [1000 40000];
vg = -2500:100:2500;


h = setup_subplots(3,1); isub = 1;
hca = h(isub); isub = isub + 1;
vdf_2D_MN.plot_plane(hca);
hca.XLabel.String = 'v_M';
hca.YLabel.String = 'v_N';

hca = h(isub); isub = isub + 1;
plotyy(hca,vdf_1D_M.depend{1}(1,:),squeeze(vdf_1D_M.data(1,:,:)),...
         vdf_2D_MN.depend{1}(1,:),sum(squeeze(vdf_2D_MN.data(1,:,:)),2))
hca.XLabel.String = 'v_M';

hca = h(isub); isub = isub + 1;
plotyy(hca,vdf_1D_N.depend{1}(1,:),squeeze(vdf_1D_N.data(1,:,:)),...
         vdf_2D_MN.depend{2}(1,:),sum(squeeze(vdf_2D_MN.data(1,:,:)),1))
hca.XLabel.String = 'v_N';

irf_plot_axis_align(h) 
if 0 % LM
h = setup_subplots(3,1); isub = 1;
hca = h(isub); isub = isub + 1;
%imagesc(hca,vdf_2D_LM.depend{1}(1,:),vdf_2D_LM.depend{1}(1,:),squeeze(vdf_2D_LM.data(1,:,:)))
vdf_2D_LM.plot_plane(hca);

hca = h(isub); isub = isub + 1;
plotyy(hca,vdf_1D_L.depend{1}(1,:),squeeze(vdf_1D_L.data(1,:,:)),...
         vdf_2D_LM.depend{1}(1,:),sum(squeeze(vdf_2D_LM.data(1,:,:)),2))
hca.XLabel.String = 'v_L';

hca = h(isub); isub = isub + 1;
plotyy(hca,vdf_1D_M.depend{1}(1,:),squeeze(vdf_1D_M.data(1,:,:)),...
         vdf_2D_LM.depend{2}(1,:),sum(squeeze(vdf_2D_LM.data(1,:,:)),1))
hca.XLabel.String = 'v_M';

irf_plot_axis_align(h)    
end
%% Plot
tint_figure = irf.tint('2017-07-11T22:33:40.00Z/2017-07-11T22:34:30.00Z');
tint_figure = tint_ifred;

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];

if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic) 
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 1 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic) 
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 0 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x,mvaVe?.y,mvaVe?.z},''comp'');',ic)
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % e psd x,y,z, 3 panels
  comp_strs = ['x','y','z'];
  comps = ['L','M','N'];
  for icomp = 1:numel(comps)
    comp = comps(icomp);
    comp_str = comp_strs(icomp);
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('f1D = ef1D%s?_elim;',comp),ic)
    irf_spectrogram(hca,f1D.specrec('velocity_1D'));  
    hca.YLim = f1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{mvaVe?.(comp_str),mvaVExB?.(comp_str).resample(mvaVe?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_e','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{e%s}',comp),'(10^3 km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
    hca.YTick = -100000:20000:100000;
    hca.YTickLabels = arrayfun(@(s) sprintf('%g',s),hca.YTick,'UniformOutput',false);
    hca.YLim = 49999*[-1 1];
  end
end

if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end

if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 1 % elim for reduced distributions
    hold(hca,'on')
    c_eval(sprintf('f1D = if1D%s?_elim;','L'),ic)
    c_eval('tselim = irf.ts_scalar(f1D.time,f1D.ancillary.energy(:,1)-f1D.ancillary.delta_energy_minus(:,1));',ic)
    c_eval('lineScpot = irf_plot(hca,tselim,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')    
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,cmap) 
end
if 1 % i psd x,y,z, 3 panels
  comp_strs = ['x','y','z'];
  comps = ['L','M','N'];
  for icomp = 1:numel(comps)
    comp = comps(icomp);
    comp_str = comp_strs(icomp);
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('f1D = if1D%s?_elim;',comp),ic)
    irf_spectrogram(hca,f1D.specrec('velocity_1D'),'lin');  
    hca.YLim = f1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{mvaVi?.(comp_str),mvaVExB?.(comp_str).resample(mvaVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(10^3 km/s)'}; 
    hca.YLabel.Interpreter = 'tex'; 
    if 1 % fitted populations
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{ts.vd});',ic)      
      hold(hca,'off')
      irf_legend(hl,{'fit1','fit2','fit3'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    %hca.YTick = -100000:20000:100000;
    %hca.YTickLabels = arrayfun(@(s) sprintf('%g',s),hca.YTick,'UniformOutput',false);
    %hca.YLim = 49999*[-1 1];
    hca.YLim = f1D.depend{1}(1,[1 end]);  
  end
end

if 1 % fit
  isub = isub + 1;
  zoomy = [zoomy isub];
      hca = irf_panel('T_fit');
      irf_plot(hca,ts.n)
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end


irf_zoom(h,'x',tint_figure)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

%% 2D distribution and pressure contributionss
time = 1;
time =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
vint = [-Inf Inf];

nshift = 12;
%tint_fit =  irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:00.00Z');
tint_fit =  irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:00.00Z');
tint_vdf =  irf_time('2017-07-11T22:34:00.00Z','utc>EpochTT');
tint_vdf = tint_fit([1 2]);

%vg = 20000:100:2500; % higher grid resolution slows down fit
%c_eval('vdf = ePDist?.tlim(time).reduce(''2D'',M,N,''vint'',vint,''scpot'',scPot?,''vg'',vg);',ic)
c_eval('dist = ePDist?.tlim(time+0.015*[-1 1]+nshift*0.03);',ic)
c_eval('scpot = scPot?.resample(dist);',ic)
vg = -70000:2000:70000;
vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
vdf = vdf([1]);

h = setup_subplots(2,1); isub = 1;

hca = h(isub); isub = isub + 1;
vdf.plot_plane(hca);
colormap(hca,pic_colors('candy4'))

hca = h(isub); isub = isub + 1;
vdf.plot_plane(hca,'off-diag-pres-cont',mvaVe1.x.resample(dist),mvaVe1.y.resample(dist));

colormap(hca,pic_colors('blue_red'))
hca.CLim = max(abs(hca.CLim))*[-1 1];
axis(hca,'square')

%% 2D distribution and pressure contributions, 2 different times

%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:04.00Z'],'utc>EpochTT')+0.06;

times = times(1):0.5:times(2);
vint = [-Inf Inf];
vg = -60000:2000:60000;
elim = [100 Inf];
h = setup_subplots(2,times.length,'vertical'); isub = 1;

for itime = 1:times.length
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h(isub); isub = isub + 1;
  vdf.plot_plane(hca);
  colormap(hca,pic_colors('candy4'))
  
  hca = h(isub); isub = isub + 1;
  vdf.plot_plane(hca,'off-diag-pres-cont',mvaVe1.x.resample(dist),mvaVe1.y.resample(dist));
  
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  %axis(hca,'square')
end

hlinks_all = linkprop(h,{'XLim','YLim'});
hlinks_f = linkprop(h(1:2:end),{'CLim'});
hlinks_p = linkprop(h(2:2:end),{'CLim'});

h(1).XLim = vg([1 end])*1e-3;
h(1).YLim = vg([1 end])*1e-3;

c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
c_eval('axis(h(?),''square'');',1:numel(h))

%% 2D distribution and pressure contributions, MN
ic = 1;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:04.00Z'],'utc>EpochTT')+0.06;

times = times(1):0.5:times(2);
vint = [-Inf Inf];
vg = -60000:2000:60000;
elim = [100 Inf];
h = setup_subplots(2,times.length,'vertical'); isub = 1;

for itime = 1:times.length
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h(isub); isub = isub + 1;
  vdf.plot_plane(hca);
  colormap(hca,pic_colors('candy4'))
  
  hca = h(isub); isub = isub + 1;
  vdf.plot_plane(hca,'off-diag-pres-cont',mvaVe1.y.resample(dist),mvaVe1.z.resample(dist));
  
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  %axis(hca,'square')
end

hlinks_all = linkprop(h,{'XLim','YLim'});
hlinks_f = linkprop(h(1:2:end),{'CLim'});
hlinks_p = linkprop(h(2:2:end),{'CLim'});

h(1).XLim = vg([1 end])*1e-3;
h(1).YLim = vg([1 end])*1e-3;

c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
c_eval('axis(h(?),''square'');',1:numel(h))

%% 2D distribution and pressure contributions, MN, X times, with locations shown
ic = 3;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
times = times + 0.25;

vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(3,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 
%delete(h1_pos);
for ip = 1:(nt*2)
  h2(ip) = subplot(3,nt,nt+ip);
end

hca = h1(1);
irf_plot(hca,mvaPe1.xy,'color','k','linewidth',1)
hca.YLabel.String = 'P_{eLM} (nPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
hmark = irf_pl_mark(hca,times.epochUnix,'k');
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;

isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',mvaVe1.x.resample(dist),mvaVe1.y.resample(dist));  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
end

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)


%% Plot data fomr simulation
% ons = PIC('/Volumes/DataRaid/Susanne-onset/data_h5/fields.h5');
% no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');
pic = ons.twcilim(65).xlim([5 25]).zlim([-2 2]);
clims = {0.01*[-1 1],0.03*[-1 1],0.01*[-1 1],0.03*[-1 1]};

pic = no02m.twpelim([16000]).xlim([95 110]).zlim([-2 2]);
clims = {0.003*[-1 1],0.3*[-1 1],0.003*[-1 1],0.3*[-1 1]};

varstrs = {'pexy','-dxPexy./ne','peyz','-dzPezy./ne'}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};
pic.plot_map(varstrs,'clim',clims,'cmap',cmaps,'A',0.2,'smooth',5)


%%
pic = no02m.twpelim([10000 24000]).xlim([95 110]).zlim(0.2*[-1 1]);
varstrs = {'vex','vey','pexy','-dxPexy./ne'}';
clims = {7*[-1 1],7*[-1 1],0.01*[-1 1],0.3*[-1 1],0.003*[-1 1],0.3*[-1 1]};
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};
pic.plot_timemap('xt',varstrs(:),'clim',clims(:),'cmap',cmaps(:),'A',0.2,'smooth',5)

% timemap



%% Plot data fomr simulation, electron ohms law for separate populations
% no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');

pic = no02m.twpelim([20000]).xlim([95 110]).zlim([-2 2]);
clims = {0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1]};

figure(51)
varstrs = {'vzdzvy([2 4 6])/100','-dzPzy([2 4 6])./n([2 4 6])','vzdzvy([4])/100','-dzPzy([4])./n([4])'}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};
pic.plot_map(varstrs,'clim',clims,'cmap',cmaps,'A',0.2,'smooth',5)

figure(52)
varstrs = {{'vzdzvy([2 4 6])/100','vzdzvy([4])/100'},{'-dzPzy([2 4 6])./n([2 4 6])','-dzPzy([4])./n([4])'}}';
ylims = {0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],0.3*[-1 1],1*[-1 1]};
varstrs = {{'vzdzvy([2])/100','vzdzvy([4])/100','vzdzvy([6])/100','(vzdzvy([2])+vzdzvy([4])+vzdzvy([6]))/100','vzdzvy([2 4 6])/100'},{'-dzPzy([2])./n([2])','-dzPzy([4])./n([4])','-dzPzy([6])./n([6])','-dzPzy([2 4 6])./n([2 4 6])'}}';
varstrs = {{'vzdzvy([2])/100','vzdzvy([4])/100','vzdzvy([6])/100'},...
           {'(vzdzvy([2])+vzdzvy([4])+vzdzvy([6]))/100','vzdzvy([2 4 6])/100'},...
           {'-dzPzy([2])./n([2])','-dzPzy([4])./n([4])','-dzPzy([6])./n([6])'},...           
           {'-dzPzy([2])./n([2])-dzPzy([4])./n([4])-dzPzy([6])./n([6])','-dzPzy([2 4 6])./n([2 4 6])'},...
           {'-(vzdzvy([2])+vzdzvy([4])+vzdzvy([6]))/100-(dzPzy([2])./n([2])+dzPzy([4])./n([4])+dzPzy([6])./n([6]))'}}';
% without hot electrons
varstrs = {{'vzdzvy([4])/100','vzdzvy([6])/100'},...
           {'(vzdzvy([4])+vzdzvy([6]))/100','vzdzvy([4 6])/100'},...
           {'-dzPzy([4])./n([4])','-dzPzy([6])./n([6])'},...           
           {'-dzPzy([4])./n([4])-dzPzy([6])./n([6])','-dzPzy([4 6])./n([4 6])'},...
           {'-(vzdzvy([4])+vzdzvy([6]))/100-(dzPzy([4])./n([4])+dzPzy([6])./n([6]))','-vzdzvy([4 6])/100-dzPzy([4 6])./n([4 6])','Ey'}}';
pic.xlim([101 101.5]+2).zlim([-1 1]).plot_line('z',varstrs,'ylim',ylims,'smooth',5)






%% 2D distribution and pressure contributions, LM, X times, comparing different satellites
ic = 1;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z';...
  '2017-07-11T22:34:03.670000000Z';...
  '2017-07-11T22:34:04.300000000Z']);
%times = times([2 3 4])+0.25;
times = times + 0.05 + 0*0.25;


vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(5,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 

%delete(h1_pos);
for ip = 1:(nt*2*2)
  h2(ip) = subplot(5,nt,nt+ip);
end

hca = h1(1);
set(hca,'colororder',mms_colors('12'))
irf_plot(hca,{mvaPe1.yz,mvaPe3.yz},'comp')
hold(hca,'on')
fhigh = 2;
irf_plot(hca,{mvaPe1.yz.filt(0,fhigh,[],5),mvaPe3.yz.filt(0,fhigh,[],5)},'comp')
hold(hca,'off')

hca.YLabel.String = 'P_{eMN} (nPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
hmark = irf_pl_mark(hca,times.epochUnix,'k');
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
%
isub = 1;
for ic = [1 3]
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca);
    hc_.Colorbar.YLabel.String = 'log_{10} f_e(v_M,v_N) (s^2/m^5)';
    colormap(hca,pic_colors('candy4'))   
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    c_eval('ve = mvaVe?;',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.y.resample(dist),ve.z.resample(dist));  
    hc_.Colorbar.YLabel.String = 'f_e(v_M,v_N)(v_M-v_M^{bulk})(v_N-v_N^{bulk}) (1/m^3)';
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
end

hlinks_all = linkprop(h2,{'XLim','YLim'});
%hlinks_f = linkprop(h2([1:nt]*[0 2]'),{'CLim'});
%hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});
hlinks_f = linkprop(h2([1:nt 2*nt+(1:nt)]),{'CLim'});
hlinks_p = linkprop(h2([(nt+1):2*nt 2*nt+((nt+1):2*nt)]),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)









