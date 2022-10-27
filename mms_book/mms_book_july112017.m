%% Load data
ic = 1;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');


% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)

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
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc



c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

%% Rotate things into new coordinate system 
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9239];
lmn = [L;M;N];


% Rotate data
c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)

c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''VExB LMN'';',ic)
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
%c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
%c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)
%mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)




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
if 0 % electron LMN, elim
  tint_efred =  irf.tint('2017-07-11T22:33:40.00Z/2017-07-11T22:34:30.00Z');
  elim = [100 40000];
  c_eval('ef1DL?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',L,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DM?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',M,''vint'',vint,''scpot'',scPot?);',ic)
  c_eval('ef1DN?_elim = ePDist?.elim(elim).tlim(tint_efred).reduce(''1D'',N,''vint'',vint,''scpot'',scPot?);',ic)
end
%% Estimating the ion components
tint_fit =  irf.tint('2017-07-11T22:33:40.00Z/2017-07-11T22:34:10.00Z');
vdf = if1DN1_elim.tlim(tint_fit);
vdf = vdf(1:5:vdf.length);
nPop = 3;
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  0,      1000e3];
[fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',0,'guessprevious',1,'X0',X0);



%% Plot
tint_figure = irf.tint('2017-07-11T22:33:40.00Z/2017-07-11T22:34:30.00Z');
tint_figure = tint_ifred;

npanels = 7;
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
    c_eval(sprintf('f1D = if1D%s?_elim;',comp),ic)
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
    irf_spectrogram(hca,f1D.specrec('velocity_1D'));  
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
