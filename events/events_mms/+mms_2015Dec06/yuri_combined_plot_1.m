%% Load data
tint = irf.tint('2015-12-06T23:38:15.00/2015-12-06T23:38:50.00Z');
ic = 1:4;
c_eval('tic; gseB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',tint); toc',ic);
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_dmpa'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; desDist? = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint); toc',ic);
c_eval('tic; disDist? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-dist'',''mms?_dis_dist_brst'',tint); toc',ic);
R  = mms.get_data('R_gse',tint);
c_eval('gseR? = TSeries(R.time,R.gseR?,''to'',1);',ic);
%%
db_info = datastore('mms_db');   
ic = 2:4;
% Electrons distribution help data
disp('Loading electron distribution help data...')
c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/12/06/mms?_fpi_brst_l2_des-dist_20151206233734_v2.1.0.cdf'']);',3:4);
c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/12/06/mms?_fpi_brst_l2_des-dist_20151206233734_v2.0.0.cdf'']);',2);
%c_eval('desDist? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_dist_brst''));',ic);
c_eval('eenergy0? = get_variable(etmpDataObj?,''mms?_des_energy0_brst'');',ic);
c_eval('eenergy1? = get_variable(etmpDataObj?,''mms?_des_energy1_brst'');',ic);
c_eval('ephi? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_phi_brst''));',ic);
c_eval('etheta? = get_variable(etmpDataObj?,''mms?_des_theta_brst'');',ic);
c_eval('estepTable? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_steptable_parity_brst''));',ic);
  
%% Set up figure
tintOverview = irf.tint('2015-12-06T23:38:22.00/2015-12-06T23:38:42.00Z');
tint = tintOverview;
%tintDist = irf.tint('2015-10-16T10:33:20.00Z',0.03);
[h1,h2] = initialize_combined_plot(10,3,2,3,'horizontal');

%% Plot timeseries 
tintOverview = irf.tint('2015-12-06T23:38:15.00/2015-12-06T23:38:50.00Z');
tint = tintOverview;

hca = irf_panel('Bx');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.99 0.9]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('abs B');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.abs.tlim(tint),dmpaB2.abs.tlim(tint),dmpaB3.abs.tlim(tint),dmpaB4.abs.tlim(tint)},'comp');
hca.YLabel.String = {'|B|','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('vx');
set(hca,'ColorOrder',mms_colors('1234'))
try; irf_plot(hca,{dbcsVe1.x.tlim(tint),dbcsVe2.x.tlim(tint),dbcsVe3.x.tlim(tint),dbcsVe4.x.tlim(tint)},'comp');
catch; irf_plot(hca,{dbcsVe2.x.tlim(tint),dbcsVe2.x.tlim(tint),dbcsVe3.x.tlim(tint),dbcsVe4.x.tlim(tint)},'comp'); end
hca.YLabel.String = {'v_{x}','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('vy');
set(hca,'ColorOrder',mms_colors('1234'))
try; irf_plot(hca,{dbcsVe1.y.tlim(tint),dbcsVe2.y.tlim(tint),dbcsVe3.y.tlim(tint),dbcsVe4.y.tlim(tint)},'comp');
catch; irf_plot(hca,{dbcsVe2.y.tlim(tint),dbcsVe2.y.tlim(tint),dbcsVe3.y.tlim(tint),dbcsVe4.y.tlim(tint)},'comp'); end
hca.YLabel.String = {'v_{y}','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('vz');
set(hca,'ColorOrder',mms_colors('1234'))
try; irf_plot(hca,{dbcsVe1.z.tlim(tint),dbcsVe2.z.tlim(tint),dbcsVe3.z.tlim(tint),dbcsVe4.z.tlim(tint)},'comp');
catch; irf_plot(hca,{dbcsVe2.z.tlim(tint),dbcsVe2.z.tlim(tint),dbcsVe3.z.tlim(tint),dbcsVe4.z.tlim(tint)},'comp'); end
hca.YLabel.String = {'v_{z}','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('abs v');
set(hca,'ColorOrder',mms_colors('1234'))
try; irf_plot(hca,{dbcsVe1.abs.tlim(tint),dbcsVe2.abs.tlim(tint),dbcsVe3.abs.tlim(tint),dbcsVe4.abs.tlim(tint)},'comp');
catch; irf_plot(hca,{dbcsVe2.abs.tlim(tint),dbcsVe2.abs.tlim(tint),dbcsVe3.abs.tlim(tint),dbcsVe4.abs.tlim(tint)},'comp'); end
hca.YLabel.String = {'|v|','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

if 1
hca = irf_panel('brst n');
set(hca,'ColorOrder',mms_colors('1234'))
try; irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
catch; irf_plot(hca,{ne2.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp'); end
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('1234'))
end

if 0
hca = irf_panel('sc Pot');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{scPot1.tlim(tint),scPot2.tlim(tint),scPot3.tlim(tint),scPot4.tlim(tint)},'comp');
hca.YLabel.String = {'sc Pot','(V)'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('1234'))
end


irf_zoom(h1(:),'x',tintOverview)
irf_zoom(h1,'y')
irf_plot_axis_align

%% Plot single time particle distributions, 4sc projections
tintDist = irf.tint('2015-12-06T23:38:31.00/2015-12-06T23:38:32.00Z');
c_eval('times = desDist?.tlim(tintDist).time;',2)


for ii = 19%:times.length;
  time = times(ii);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  vlim = 15*1e3; % velocity limit for projection plots
  elevlim = 20; % elevation limit to be included in projection
  strCMap = 'jet'; % colormap  
  projclim = [-1 4.5];
  palim = [1e-3 1e6]; % pitch angle distribution plot y-limit
  skymapEnergy = [65 278]; % skymap energies
   
  c_eval('dist = desDist?;',ic)
   
  isub = 1; 
  hca = h2(isub); isub = isub + 1;
  plot_lmn3D(hca,gseR1,gseR2,gseR3,gseR4,[1 0 0;0 1 0;0 0 1],{'X','Y','Z'})
  view(hca,[0 -1 0])
  axis(hca,'square')

  hca = h2(isub); isub = isub + 1;
  plot_lmn3D(hca,gseR1,gseR2,gseR3,gseR4,[1 0 0;0 1 0;0 0 1],{'X','Y','Z'})  
  view(hca,[0 0 1])
  axis(hca,'square')
  
  for ic = 2:4 % mms1 data missing
    hca = h2(isub); isub = isub + 1;
    c_eval('Ve0 = dbcsVe?.resample(time).data;',ic); 
    hatVe0 = double(irf_norm(Ve0));    

    % Get mean magnetic field direction
    tint = time + 0.5*0.03*[-1 1]; % averge over electron distribution
    c_eval('B0 = mean(dmpaB?.resample(dslE?.tlim(tint).time).data);',ic); 
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?.tlim(tint).data);',ic); 
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);

    vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e'};%0;hatVe0,'V_e'};

    % Projection coordinate system
    x = hatB0;
    y = hatExB0;
    z = cross(x,y);
    vlabel = {'v_{B,dir}','v_{ExB,dir}','v_{E,dir}'};
    
    % Plot project ion onto a plane
    c_eval('phi = ephi?; theta = etheta?; stepTable = estepTable?; energy0 = eenergy0?; energy1 = eenergy1?;',ic)      
    mms.plot_projection(hca,dist,phi,theta,stepTable,energy0,energy1,...
                'tint',times(ii),'xyz',[x;y;z],'elevationlim',elevlim,...
                'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',vlabel);
    hca.Title.String = irf_ssub('MMS ?',ic);
    colormap(hca,strCMap)
  end
    
  pause(0.1)
  %cn.print([irf_ssub('BvnP_psds_mms?_',ic) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

%% Plot single time particle distributions, 1sc projections, pitch angles, skymap
tintDist = irf.tint('2015-12-06T23:38:31.00/2015-12-06T23:38:32.00Z');
c_eval('times = desDist?.tlim(tintDist).time;',ic)
ic  = 4;
%times = irf.tint('2015-12-06T23:38:31.60/2015-12-06T23:38:32.00Z');
for ii = 19%:times.length;
  time = times(ii);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  vlim = 15*1e3; % velocity limit for projection plots
  elevlim = 20; % elevation limit to be included in projection
  strCMap = 'jet'; % colormap  
  projclim = [-1 4.5];
  palim = [1e-3 1e6]; % pitch angle distribution plot y-limit
  skymapEnergy = [65 278]; % skymap energies
   
  c_eval('dist = desDist?;',ic)
 
  c_eval('Ve0 = dbcsVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
   
  % Get mean magnetic field direction
  tint = time + 0.5*0.03*[-1 1]; % averge over electron distribution
  c_eval('B0 = mean(gseB?.resample(dslE?.tlim(tint).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = mean(dslE?.tlim(tint).data);',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
  
  vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e'};%0;hatVe0,'V_e'};

  % Projection coordinate system
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  
  isub = 1;
  
  % Plot psd 0 90 180
  hca = h2(isub); isub = isub + 1;
  psdtint = times(ii);
  c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?,''tint'',psdtint,''scPot'',scPot?,''ylim'',palim,''energies'',skymapEnergy);',ic)
  TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
  %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
  %hca.Title.String = irf_ssub('C?',ic);
  hca.Title.String = TitleStr;
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',times(ii),''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',times(ii),''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};

  % Plot project ion onto a plane
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)

  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
    
  pause(0.1)
  %cn.print([irf_ssub('BvnP_psds_mms?_',ic) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end