% Scattering and trapping figure
% Load data
tint = irf.tint('2015-10-16T10:30:00.00Z/2015-10-16T10:30:40.00Z');
tint = irf.tint('2015-10-16T10:32:40.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere

tint_figure = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z');

mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');

ic = 1;
%%
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)

c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('ePitch? = ePDist?.tlim(tint_figure).pitchangles(dmpaB?);',ic)

%% Magnetic curvature
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');


%% Figure

npanels = 3;
nrows = 2;
ncols = 1;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.5,'vertical'); % horizontal


elim_pitch = [190 200];
% TSeries plots
if 1 % B   
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B^{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 0 % E   
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E^{GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % vExB   
  hca = irf_panel('vExB');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)
  hca.YLabel.String = {'v_{ExB}^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 1 % iPitch  
  hca = irf_panel('ePitch');
  c_eval('iPitch = ePitch?.elim(elim_pitch);',ic)
  irf_spectrogram(hca,iPitch.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,pic_colors('pasteljet'))
end


irf_plot_axis_align
irf_zoom(h1,'x',tint_figure)
irf_zoom(h1(1),'y')

%%
isub = 1;
if exist('hmark'); delete(hmark); end
times_dist =  EpochTT(['2015-10-16T10:33:30.326Z']) + 0*0.03*[1];
%times_dist =  EpochTT(['2015-10-16T10:33:30.386Z']) + 0*0.03*[1];

c_eval('scpot = scPot?;',ic)
vint = 4e3*[-1 1];
vlim = 20;
for id = 1:times_dist.length  
  
  dt = 0.03*0.5;
  c_eval('v1 = mean(dmpaB?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  c_eval('v2 = mean(dslE?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  v2 = cross(v1,cross(v2,v1)); v2 = v2/norm(v2); % E
  v3 = cross(v1,v2); % -ExB
  vg = (-40:1:40)*1e3;
  %vlim = 70;
  if 1 % vB,vExB
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',-v3,v1,'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.YLabel.String = 'v_B (10^3 km/s)';
    hca.XLabel.String = 'v_{ExB} (10^3 km/s)';
    colormap(hca,pic_colors('candy4'))
  end
  if 0 % vE,vExB
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v2,-v3,'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_E (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  colormap(hca,pic_colors('candy4'))
  end
  if 0 % vB,vE
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{E} (10^3 km/s)';
    colormap(hca,pic_colors('candy4'))
  end
  if 0 % x,y
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',[1 0 0],[0 1 0],'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_x (10^3 km/s)';
    hca.YLabel.String = 'v_{y} (10^3 km/s)';
    colormap(hca,pic_colors('candy4'))
  end
  if 0 % y,z
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',[0 1 0],[0 0 1],'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_y (10^3 km/s)';
    hca.YLabel.String = 'v_{z} (10^3 km/s)';
  colormap(hca,pic_colors('candy4'))
  end
  if 0 % z,x
    hca = h2(isub); isub = isub + 1;
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',[0 0 1],[1 0 0],'vg',vg,'vint',vint,'scpot',scpot.resample(ddist));

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_z (10^3 km/s)';
    hca.YLabel.String = 'v_{x} (10^3 km/s)';
  end
  colormap(hca,pic_colors('candy4'))
end
for id = 1:times_dist.length
  hca = h2(isub); isub = isub + 1;
  dt = 0.05;
  
  c_eval('ddist = ePitch?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
  ddist = ddist(1);
  times_dist_exact_eu(id) = ddist.time.epochUnix;

  hh = ddist.plot_pad_polar(hca,'scpot',scpot.resample(ddist).data);
  %axis(hca,'square')
%   hca.XLim = 0.99*vlim*[-1 1];
%   hca.YLim = 0.99*vlim*[-1 1];
  hca.Layer = 'top';
  %hca.XLabel.String = 'v_E (10^3 km/s)';
  %hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  
  hca.XLim = [-3.6 0];
  hca.YLim = [-3.6 3.6];
  shading(hca,'flat')
  colormap(hca,pic_colors('candy4'))
end

%c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
hmark = irf_pl_mark(h1,tocolumn(times_dist_exact_eu),'r');
c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))
