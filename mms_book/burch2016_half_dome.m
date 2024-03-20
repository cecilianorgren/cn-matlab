%% Load data
% units = irf_units;
irf.log('critical')
ic = 2;

localuser = 'cno062';
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
%mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
mms.db_init('local_file_db',['/Volumes/mms']);
db_info = datastore('mms_db');   

tint_all = irf.tint('2015-10-16T13:05:24.00Z/2015-10-16T13:07:30.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);
iFile = 1;
tint = [files(iFile).start files(iFile).stop] + 0*[1 -1];

tic
% Magnetic field
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

% Electric field
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);

% Spacecraft potential
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);

% Velocity
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)

% Spacecraft position
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic)

% Distributions
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data

toc

%% Plot distributions
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');
%times = times + 0.090;
times = times(2:end);
nt = times.length;

vint_red = [-2000 2000];
vint = [3500 5000];
vint = [4000 5000] - 500;
vint = [3000 6000];
eint = [50 90];
vlim_red = [-10 10];
eint = units.me*(vint*1e3).^2/2/units.eV;
c_eval('PDist = ePDist?;',ic)
c_eval('dmpaB = dmpaB?;',ic)
c_eval('dslE = dslE?;',ic)
c_eval('scPot = scPot?;',ic)

nrows = 3;
ncols = 4;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

for it = 1:nt % perp perp
  hca = h(isub); isub = isub + 1;
  time = times(it);
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  %scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);
  fred = pdist.reduce('2D',ExB,Eperp,'scpot',scPot,'vint',vint_red,'vg',-30000:500:30000);
  fred.plot_plane(hca);
  %colormap(hca,irf_colormap('waterfall'))
  hca.XLim = vlim_red;
  hca.YLim = vlim_red;
  hca.YLabel.String = 'v_{E_\perp}';
  hca.XLabel.String = 'v_{ExB}';
  if 1 % plot vint
    hold(hca,'on')
    plot(hca,1e-3*vint(1)*cosd(0:360),1e-3*vint(1)*sind(0:360),'k')
    plot(hca,1e-3*vint(2)*cosd(0:360),1e-3*vint(2)*sind(0:360),'k')
    hold(hca,'off')    
  end
  if 1 % plot vint_red
    hold(hca,'on')
    plot(hca,1e-3*vint_red(1)*[1 1],[-30 30],'k')
    plot(hca,1e-3*vint_red(2)*[1 1],[-30 30],'k')
    hold(hca,'off')    
  end
  axis(hca,'square')
end

for it = 1:nt % par perp
  hca = h(isub); isub = isub + 1;
  time = times(it);
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  %scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);
  fred = pdist.reduce('2D',B,ExB,'scpot',scPot,'vint',vint_red,'vg',-30000:500:30000);
  fred.plot_plane(hca);
  %colormap(hca,irf_colormap('waterfall'))
  hca.XLim = vlim_red;
  hca.YLim = vlim_red;
  hca.XLabel.String = 'v_{B}';
  hca.YLabel.String = 'v_{ExB}';
  if 1 % plot vint_skymap
    hold(hca,'on')
    plot(hca,1e-3*vint(1)*cosd(0:360),1e-3*vint(1)*sind(0:360),'k')
    plot(hca,1e-3*vint(2)*cosd(0:360),1e-3*vint(2)*sind(0:360),'k')
    hold(hca,'off')    
  end
  if 1 % plot vint_red
    hold(hca,'on')
    plot(hca,1e-3*vint_red(1)*[1 1],[-30 30],'k')
    plot(hca,1e-3*vint_red(2)*[1 1],[-30 30],'k')
    hold(hca,'off')    
  end
  axis(hca,'square')
end

for it = []%1:nt % spherical skymap
  hca = h(isub); isub = isub + 1;
  time = times(it);
  pdist = PDist.elim(eint).tlim(time+0.499*0.03*[-1 1]);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);
  [ax,hcb] = mms.plot_skymap(hca,pdist,'vectors',{B,'B';E,'E';Eperp,'E_\perp';ExB,'ExB'},'smooth',1);
  colormap(hca,irf_colormap('waterfall'))
end

for it = 1:nt % flat skymap
  %%  
  hca = h(isub); isub = isub + 1;
  time = times(it);
  pdist = PDist.elim(eint).tlim(time+0.499*0.03*[-1 1]);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);
  [ax,hcb] = mms.plot_skymap(hca,pdist,'flat','vectors',{B,'B';E,'E';Eperp,'E_\perp';ExB,'ExB'},'smooth',2,'pitchangle',B,[45 90 135]);
  colormap(hca,irf_colormap('waterfall'))
end

if 1
  hlinks1 = linkprop(h(1:2*ncols),{'CLim'});
  h(1).CLim = [-9 -7.];
  hlinks2 = linkprop(h(2*ncols+1:end),{'CLim'});
  hlinks3 = linkprop(h(2*ncols+(1:ncols)),{'View'});
else
  hlinks1 = linkprop(h(1:2*nrows),{'CLim'});
  h(1).CLim = [-9 -7.];
  hlinks2 = linkprop(h(2*nrows+1:end),{'CLim'});
  hlinks3 = linkprop(h(2*nrows+(1:nrows)),{'View'});

end

%c_eval('h(?).XLim = 15*[-1 1]; h(?).YLim = 15*[-1 1];',1:nrows)
%c_eval('h(?).XLim = 15*[-1 1]; h(?).YLim = 15*[-1 1];',1:nrows)
h(end).CLim = [0 2]*1e-26;
drawnow
colormap([0.8 0.8 0.8; pic_colors('candy4')])





%% Plot isosurfaces of distribution
units = irf_units;
ic = 4;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');
times = times + 0.090;
nt = times.length;

vint_red = [-1000 1000];
vint = [3500 5000];
vint = [3000 5000] + 500;
eint = [50 90];
vlim_red = [-10 10];
eint = units.me*(vint*1e3).^2/2/units.eV;
c_eval('PDist = ePDist?;',ic)
c_eval('dmpaB = dmpaB?;',ic)
c_eval('dslE = dslE?;',ic)
c_eval('scPot = scPot?;',ic)


nrows = 5;
ncols = 1;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

for it = 1:nt % par perp
  hca = h(isub); isub = isub + 1;
  time = times(it);
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  %scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);

  F = squeeze(pdist.data);
  Flev = prctile(F(:),95);
  [faces,verts] = extractIsosurface(F,Flev);

  figure
  p = patch(Faces=faces,Vertices=verts);
  isonormals(V,p)
  view(3)
  set(p,FaceColor=[0.5 1 0.5])
  set(p,EdgeColor="none")
  camlight
  lighting gouraud

  colormap(hca,irf_colormap('waterfall'))
  %hca.XLim = vlim_red;
  %hca.YLim = vlim_red;
  %hca.XLabel.String = 'v_{B}';
  %hca.YLabel.String = 'v_{ExB}';

  axis(hca,'square')
end











