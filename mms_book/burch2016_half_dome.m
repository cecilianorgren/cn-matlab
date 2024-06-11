%% Load data
% units = irf_units;
irf.log('critical')
ic = 2:3:4;

localuser = 'cno062';
%localuser = 'cecilia';
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
%mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
%mms.db_init('local_file_db',['/Volumes/mms']);
db_info = datastore('mms_db');   
%
tint_all = irf.tint('2015-10-16T13:00:24.00Z/2015-10-16T13:09:30.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);
iFile = 1;
tint = [files(iFile).start files(iFile).stop] + 0*[1 -1];

tic
% Magnetic field
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

% Electric field
c_eval('gseE?fast = mms.get_data(''E_gse_edp_fast_l2'',tint,?);',ic);
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

Lgse = [0.3665, -0.1201, 0.9226];
Mgse = [0.5694, -0.7553, -0.3245];
Ngse = [0.7358, 0.6443, -0.2084];
Tgse = [Lgse; Mgse; Ngse];
c_eval('tsLgse? = irf.ts_vec_xyz(ePDist?.time,repmat(Lgse,ePDist?.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(ePDist?.time,repmat(Mgse,ePDist?.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(ePDist?.time,repmat(Ngse,ePDist?.length,1));',ic)

c_eval('defatt = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)
c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt,-1);',ic)

%% Plot distributions
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');
%times = times + 0.090;
times = times(1:end);
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
ncols = 5;
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
  [ax,hcb] = mms.plot_skymap(hca,pdist,'flat','vectors',{B,'B';E,'E';Eperp,'E_\perp';ExB,'ExB'},'smooth',1,'pitchangle',B,[45 90 135]);
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
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');
%times = times + 0.090;
times = times(1:3);
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


%if exist('h','var'); delete(h)

nrows = 1;
ncols = nt;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;
%nt = 1;

set(gcf,'defaultLegendAutoUpdate','off');

for it = 1:nt % 3D isosurface
    
  hca = h(isub); isub = isub + 1;
  time = times(it);
  elimit = 40;
  elimit = 35;
  %pdist = PDist.elim([elimit Inf]).tlim(time+0.499*0.03*[-1 1]);  
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  pdist.data(:,1:5,:,:) = 3e-27;


  scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); B = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); E = E/norm(E);
  Eperp = cross(B,cross(E,B));
  ExB = cross(E,B);
  
  nSmooth = 3;
  if 1
    iso_values = [2e-27, 6e-27];
    set(hca,'colororder',[1 0.5 0.5; 0.5 1 0.5]) 
    T = [L;M;N];
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',T);
    
    hs.Patch(1).FaceAlpha = 0.3;
    hs.Patch(2).FaceAlpha = 0.9;
    
    legs = arrayfun(@(x)sprintf('%g',x),iso_values,'UniformOutput',false);
    hleg = legend(hs.Patch,legs,'location','northeast');
    hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
  else
    F = squeeze(pdist.data);
    [VX,VY,VZ] = pdist.v('scpot',scP*0);
    VX = squeeze(VX);
    VY = squeeze(VY);
    VZ = squeeze(VZ);
    if 1 % make circular in azimuth
      %F = cat(2,[F(:,end,:),F,F(:,1,:)]);
      %VX = cat(2,[VX(:,end,:),VX,VX(:,1,:)]);
      %VY = cat(2,[VY(:,end,:),VY,VY(:,1,:)]);
      %VZ = cat(2,[VZ(:,end,:),VZ,VZ(:,1,:)]);
      F = cat(2,[F,F(:,1,:)]);
      VX = cat(2,[VX,VX(:,1,:)]);
      VY = cat(2,[VY,VY(:,1,:)]);
      VZ = cat(2,[VZ,VZ(:,1,:)]);
    end
    if 0 % interpolate to close circle
      
    end    
    if 0 % Rebin to square grid, for smoothing purposes
      MP = pdist.macroparticles('skipzero',1,'ntot',1e6);
      vmax = 100e3;
      dv = 1000;
      vx = -vmax:dv:vmax;
      vy = -vmax:dv:vmax;
      vz = -vmax:dv:vmax;
      [VXg,VYg,VZg] = ndgrid(vx,vy,vz);
      %[count edges mid loc] = histcn([VX(:),VY(:),VZ(:)],vx,vy,vz,'fun',@sum,'accumdata',MP.dn);  
      [count edges mid loc] = histcn([MP.vx,MP.vy,MP.vz],vx,vy,vz,'fun',@sum,'accumdata',MP.dn);  
    
      F = count;
      VX = mid{1};
      VY = mid{2};
      VZ = mid{3};
    end
  %    
    F = smooth3(F,'box',nSmooth);
    Flev = prctile(F(:),85);
    Flev = 1e-26;
    
    Flev = 1e-27;
    Flev = prctile(F(:),85);
    Flev = 6e-27;
    s = isosurface(VX,VY,VZ,F,Flev);
    %s = isocaps(VX,VY,VZ,F,Flev);
    p = patch(hca,Faces=s.faces,Vertices=s.vertices);
    set(p,FaceColor=[0.5 1 0.5])
    set(p,EdgeColor="none")
    set(p,FaceAlpha=0.8)
    camlight(hca,0,-45)
    lighting(hca,'gouraud')

    if 1
      hold(hca,'on')
      Flev2 = 1e-27;
      s2 = isosurface(VX,VY,VZ,F,Flev2);
      p2 = patch(hca,Faces=s2.faces,Vertices=s2.vertices);
      set(p2,FaceColor=[1 0.5 0.5])
      set(p2,EdgeColor="none")
      set(p2,FaceAlpha=0.3)  
      hold(hca,'off')
    end
    set(hca,'colororder',[0.5 1 0.5; 1 0.5 0.5])
    irf_legend(hca,{sprintf('f = %g',Flev),sprintf('f = %g',Flev2)}',[0.98 0.98],'backgroundcolor',[1 1 1],'fontsize',12)

    if 0
      pdist_low = PDist.elim([0 elimit]).tlim(time+0.499*0.03*[-1 1]);  
      F = squeeze(pdist_low.data);
      F = smooth3(F,'box',nSmooth);
      [VX,VY,VZ] = pdist_low.v;
      VX = squeeze(VX);
      VY = squeeze(VY);
      VZ = squeeze(VZ);
      hold(hca,'on')
      s2 = isosurface(VX,VY,VZ,F,1e-27);
      p2 = patch(hca,Faces=s2.faces,Vertices=s2.vertices);
      set(p2,FaceColor=[0.5 0.5 1])
      set(p2,EdgeColor="none")
      set(p2,FaceAlpha=1)  
      hold(hca,'off')
    end
  
    %[faces,verts] = extractIsosurface(F,Flev);
  

    view(hca,[1 1 1])
    %%
    % figure
    % p = patch(Faces=faces,Vertices=verts);
    % isonormals(V,p)
    % view(3)
    % set(p,FaceColor=[0.5 1 0.5])
    % set(p,EdgeColor="none")
    % camlight
    % lighting gouraud
  
    %colormap(hca,irf_colormap('waterfall'))
    %hca.XLim = vlim_red;
    %hca.YLim = vlim_red;
    hca.XLabel.String = 'v_{x} (km/s)';
    hca.YLabel.String = 'v_{y} (km/s)';
    hca.ZLabel.String = 'v_{z} (km/s)';  
  
    axis(hca,'square')
    axis(hca,'equal')
    grid(hca,'on')
    hca.Title.String = {'Box averaged over',sprintf('%g (spherical) datapoints',nSmooth),pdist.time.utc};
  end
  
  if 1 % plot E, ExB, B
    %%
    hold(hca,'on')
    qmax = 10000; 
    nE = E/norm(E);
    nEperp = Eperp/norm(Eperp);
    nB = B/norm(B);
    nExB = ExB/norm(ExB);
    colors_quiv_B = [0 0 0];
    colors_quiv_E = [1 0 0];
    colors_quiv_Eperp = [1 0.5 0.5];
    colors_quiv_ExB = [0 0.5 1]; 

    if 1
      nE = nE*hs.Rotation';
      nExB = nExB*hs.Rotation';
      nB = nB*hs.Rotation';
      nEperp = nEperp*hs.Rotation';
    end
    
    fontsize = 13;
    fontweight = 'bold';
    quiver3(hca,qmax*nB(1)*-1,qmax*nB(2)*-1,qmax*nB(3)*-1,qmax*nB(1)*2,qmax*nB(2)*2,qmax*nB(3)*2,'color',colors_quiv_B,'linewidth',2)
    text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,'B','color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*nE(1)*-1,qmax*nE(2)*-1,qmax*nE(3)*-1,qmax*nE(1)*2,qmax*nE(2)*2,qmax*nE(3)*2,'color',colors_quiv_E,'linewidth',2)
    text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,'E','color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nEperp(1)*-1,qmax*nEperp(2)*-1,qmax*nEperp(3)*-1,qmax*nEperp(1)*2,qmax*nEperp(2)*2,qmax*nEperp(3)*2,'color',colors_quiv_Eperp,'linewidth',2)
    text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,'E_\perp','color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*nExB(1)*-1,qmax*nExB(2)*-1,qmax*nExB(3)*-1,qmax*nExB(1)*2,qmax*nExB(2)*2,qmax*nExB(3)*2,'color',colors_quiv_ExB,'linewidth',2)    
    text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,'E\times{B}','color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')
    %set(hca,'colororder',colors_quiv)
    %irf_legend(hca,{'B','E','E_\perp','E\times{B}'},[0.98 0.98])
  end
  if 1 % plot LMN (obs in GSE!!)
    %%
     if 1
      L_ = L*hs.Rotation';
      M_ = M*hs.Rotation';
      N_ = N*hs.Rotation';
      
     end

    hold(hca,'on')
    colors_quiv_LMN = [0.5 0.5 0.5];
    quiver3(hca,qmax*L_(1)*-1,qmax*L_(2)*-1,qmax*L_(3)*-1,qmax*L_(1)*2,qmax*L_(2)*2,qmax*L_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*L_(1)*1,qmax*L_(2)*1,qmax*L_(3)*1,'L','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*M_(1)*-1,qmax*M_(2)*-1,qmax*M_(3)*-1,qmax*M_(1)*2,qmax*M_(2)*2,qmax*M_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*M_(1)*1,qmax*M_(2)*1,qmax*M_(3)*1,'M','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*N_(1)*-1,qmax*N_(2)*-1,qmax*N_(3)*-1,qmax*N_(1)*2,qmax*N_(2)*2,qmax*N_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*N_(1)*1,qmax*N_(2)*1,qmax*N_(3)*1,'N','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')

  end
end


hlinks = linkprop(h,{'view'});
c_eval('h(?).FontSize = 12;',1:numel(h))

view([0 0 1])

if 1
  %%
  hlight = findobj(h,'type','light');
  delete(hlight)
  for ip = 1:nt
    %camlight(h(ip),45,45);
    camlight(h(ip),45,10);
    camlight(h(ip),45,190);
  end
end

if 0
  %%
  c_eval('h(?).XTick = [-8:2:8]*1e3;',1:numel(h))
  c_eval('h(?).YTick = [-8:2:8]*1e3;',1:numel(h))
  c_eval('h(?).ZTick = [-8:2:8]*1e3;',1:numel(h))
  hold(gca,'on')
  % DRAW AXIS LINEs
  plot3(get(gca,'XLim'),[0 0],[0 0],'k');
  plot3([0 0],[0 0],get(gca,'ZLim'),'k');
  plot3([0 0],get(gca,'YLim'),[0 0],'k');
  % GET TICKS
  X=get(gca,'Xtick');
  Y=get(gca,'Ytick');
  Z=get(gca,'Ztick');
  % GET LABELS
  XL=get(gca,'XtickLabel');
  YL=get(gca,'YtickLabel');
  ZL=get(gca,'ZtickLabel');
  % REMOVE TICKS
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  set(gca,'Ztick',[]);
  % GET OFFSETS
  Xoff=diff(get(gca,'XLim'))./30;
  Yoff=diff(get(gca,'YLim'))./30;
  Zoff=diff(get(gca,'ZLim'))./30;
  % DRAW TICKS
  %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
  for i=1:length(X)
     plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
  end;
  for i=1:length(Y)
     plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
  end;
  for i=1:length(Z)
     plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
  end;
  % DRAW LABELS
  text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
  text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
  text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);
  hold(gca,'off')
end
%% Plot isosurfaces of distribution, PDist, cleaned, also with TSeries of locations
units = irf_units;
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z';...
             '2015-10-16T13:07:02.310Z'];
times = irf_time(times_utc,'utc>EpochTT');

%times = EpochTT(['2017-07-11T22:34:01.300000000Z';
%  '2017-07-11T22:34:01.800000000Z';...
%  '2017-07-11T22:34:02.300000000Z';...
%  '2017-07-11T22:34:02.800000000Z';...
%  '2017-07-11T22:34:03.300000000Z']);
%times = times + 0.090;
times = times(1:3);

tint_plot = [times.start times.stop] + [-0.5 0.5];
nt = times.length;

vint_red = [-1000 1000];
vint = [3500 5000];
vint = [3000 5000] + 500;
eint = [50 90];
eint = [50 90];
vlim_red = [-10 10];
eint = units.me*(vint*1e3).^2/2/units.eV;
c_eval('PDist = ePDist?;',ic)
c_eval('dmpaB = dmpaB?;',ic)
c_eval('dslE = dslE?;',ic)
c_eval('gseE = gseE?;',ic)
c_eval('scPot = scPot?;',ic)
c_eval('gseVe = gseVe?;',ic)


%if exist('h','var'); delete(h)
fontsize = 16;

nrows = 1;
ncols = nt;
%h = setup_subplots(nrows,ncols,'vertical');
%[h1,h] = initialize_combined_plot('topbottom',2,nrows,ncols,0.5,'horizontal');
[h1,h] = initialize_combined_plot('topbottom',2,nrows,ncols,0.5,'horizontal');

isub = 1;

hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('plB = gseB?*Tgse'';',ic)
irf_plot(hca,plB)
hca.YLabel.String = 'B (nT)';  hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'L ','M ','N '},[0.98 0.98])

hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyz'))
c_eval('plVe = gseVe?*Tgse'';',ic)
irf_plot(hca,plVe)
set(hca,'ColorOrder',mms_colors('xyz'))
hold(hca,'on')
%irf_plot(hca,plE.resample(PDist,'mean'))
hold(hca,'off')
hca.YLabel.String = 'v_e (km/s)';  hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L ','M ','N '},[0.98 0.98])

if 0
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('plE = gseE?*Tgse'';',ic)
  irf_plot(hca,plE)
  set(hca,'ColorOrder',mms_colors('xyz'))
  hold(hca,'on')
  %irf_plot(hca,plE.resample(PDist,'mean'))
  hold(hca,'off')
  hca.YLabel.String = 'E (mV/m)';  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L ','M ','N '},[0.98 0.98])
end

if 0 % electric field
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('plE = gseE?fast*Tgse'';',ic)
  c_eval('plE_alt = gseE?*Tgse'';',ic)
  %plE = plE.resample(PDist.time+0.015);
  plE = plE.resample(PDist.time);
  plE_alt = plE_alt.resample(plE,'mean');
  irf_plot(hca,plE)
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,plE_alt)
  hold(hca,'off')
  hca.YLabel.String = 'E_{LMN}^{fast} (mV/m)';  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98])
end
if 0
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('plVe = gseVe?*Tgse'';',ic)
  irf_plot(hca,plVe)
  hca.YLabel.String = 'V_{e,LMN} (mV/m)';  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98])
end

h1(end).XTickLabelRotation = 0;

irf_zoom(h1,'x',tint_plot)
irf_zoom(h1,'y')

colors = get(gca,'colororder');
colors = [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.5];
for it = 1:nt % mark locations of distributions
  for ip = 1:numel(h1)
    c_eval('tcenter = ePDist?.tlim(times(it)+0.03*0.49*[-1 1]).time;',ic)
    time_int = tcenter+0.03*0.5*[-1 1];
    hmark = irf_pl_mark(h1(ip),time_int',colors(it,:),'facealpha',0.5);
  end
end


isub = 1;
%nt = 1;

set(gcf,'defaultLegendAutoUpdate','off');

for it = 1:nt % 3D isosurface
    
  hca = h(isub); isub = isub + 1;
  time = times(it);
  time = time;
  elimit = 40;
  elimit = 35;
  elimit = 30;
  pdist = PDist.elim([elimit Inf]).tlim(time+0.499*0.03*[-1 1]);  
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  pdist.data(:,1:5,:,:) = 3e-27;

  tint_tmp = time+0.5*0.03*[-1 1];
  scP = scPot.tlim(tint_tmp); scP = mean(scP.data,1);
  B = dmpaB.tlim(tint_tmp); B = mean(B.data,1); Bn = B/norm(B);
  E = dslE.tlim(tint_tmp); E = mean(E.data,1); En = E/norm(E);
  %E = gseE.tlim(tint_tmp); E = mean(E.data,1); En = E/norm(E);
  Eperp = cross(Bn,cross(E,Bn)); Eperpn = Eperp\norm(Eperp);
  ExB = cross(E,B); ExBn = Eperp\norm(Eperp);
  
  nSmooth = 3;
  
  iso_values = [0.5e-27, 6e-27];
  iso_values = [0.5e-27, 6e-27];
  %iso_values = [1.5e-30 22e-30];
  iso_values = iso_values(2); 
  set(hca,'colororder',[1 0.5 0.5; 0.5 1 0.5]) 
  set(hca,'colororder',[0.2 0.5 1; 1 0.5 0.5; 0.5 1 0.5]) 
  %set(hca,'colororder',colors(it,:)) 
  %Tgse = [Lgse;Mgse;Ngse];
  c_eval('Tdsl = [tsLdsl?.resample(pdist).data;tsMdsl?.resample(pdist).data;tsNdsl?.resample(pdist).data];',ic)
  hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);
  %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',Tdsl);
  
  hs.Patch(1).FaceAlpha = 0.9;  
  %hs.Patch(1).FaceAlpha = 0.3;
  %hs.Patch(2).FaceAlpha = 0.9;
  
  legs = arrayfun(@(x)sprintf('%g',x),hs.Values,'UniformOutput',false);
  hleg = legend(hs.Patch,legs,'location','northeast');
  hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
  
  
  if 1 % plot E, ExB, B    
    hold(hca,'on')
    qmax = 10000; 
    nE = E/norm(E);
    nEperp = Eperp/norm(Eperp);
    nB = B/norm(B);
    nExB = ExB/norm(ExB);
    colors_quiv_B = [0 0 0];
    colors_quiv_E = [1 0 0];
    colors_quiv_Eperp = [1 0.5 0.5];
    colors_quiv_ExB = [0 0.5 1]; 

    if 1
      nE = nE*hs.Rotation';
      nExB = nExB*hs.Rotation';
      nB = nB*hs.Rotation';
      nEperp = nEperp*hs.Rotation';
    end
    
    fontsize = 13;
    fontweight = 'bold';
    quiver3(hca,qmax*nB(1)*-1,qmax*nB(2)*-1,qmax*nB(3)*-1,qmax*nB(1)*2,qmax*nB(2)*2,qmax*nB(3)*2,'color',colors_quiv_B,'linewidth',2)
    %text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,'B','color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,sprintf('B=%.0fnT',norm(B)),'color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*nE(1)*-1,qmax*nE(2)*-1,qmax*nE(3)*-1,qmax*nE(1)*2,qmax*nE(2)*2,qmax*nE(3)*2,'color',colors_quiv_E,'linewidth',2)
    %text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,'E','color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,sprintf('E=%.0fmV/m',norm(E)),'color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nEperp(1)*-1,qmax*nEperp(2)*-1,qmax*nEperp(3)*-1,qmax*nEperp(1)*2,qmax*nEperp(2)*2,qmax*nEperp(3)*2,'color',colors_quiv_Eperp,'linewidth',2)
    %text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,'E_\perp','color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,['E_\perp' sprintf('=%.0fmV/m',norm(Eperp))],'color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nExB(1)*-1,qmax*nExB(2)*-1,qmax*nExB(3)*-1,qmax*nExB(1)*2,qmax*nExB(2)*2,qmax*nExB(3)*2,'color',colors_quiv_ExB,'linewidth',2)    
    text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,'E\times{B}','color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)
    %text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,sprintf('E\times{B}=%.0fmV/m',norm(Eperp)),'color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')
    %set(hca,'colororder',colors_quiv)
    %irf_legend(hca,{'B','E','E_\perp','E\times{B}'},[0.98 0.98])
  end
  if 0 % plot LMN (obs in GSE!!)
    %%
     if 1
      L_ = L*hs.Rotation';
      M_ = M*hs.Rotation';
      N_ = N*hs.Rotation';
      
     end

    hold(hca,'on')
    colors_quiv_LMN = [0.5 0.5 0.5];
    quiver3(hca,qmax*L_(1)*-1,qmax*L_(2)*-1,qmax*L_(3)*-1,qmax*L_(1)*2,qmax*L_(2)*2,qmax*L_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*L_(1)*1,qmax*L_(2)*1,qmax*L_(3)*1,'L','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*M_(1)*-1,qmax*M_(2)*-1,qmax*M_(3)*-1,qmax*M_(1)*2,qmax*M_(2)*2,qmax*M_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*M_(1)*1,qmax*M_(2)*1,qmax*M_(3)*1,'M','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*N_(1)*-1,qmax*N_(2)*-1,qmax*N_(3)*-1,qmax*N_(1)*2,qmax*N_(2)*2,qmax*N_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*N_(1)*1,qmax*N_(2)*1,qmax*N_(3)*1,'N','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')

  end
end

for ip = 1:numel(h)
  h(ip).XLabel.String = 'v_L (km/s)';
  h(ip).YLabel.String = 'v_M (km/s)';
  h(ip).ZLabel.String = 'v_N (km/s)';
  vlim = [-9 9]*1e3;
  vlim = 7.3*[-1 1]*1e3;
  h(ip).XLim = vlim;
  h(ip).YLim = vlim;
  h(ip).ZLim = vlim;
  h(ip).XTick = [-10:5:10]*1e3;
  h(ip).YTick = [-10:5:10]*1e3;
  h(ip).ZTick = [-10:5:10]*1e3;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).ZGrid = 'on';
  h(ip).Box = 'on';
end

hlinks = linkprop(h,{'view'});
c_eval('h(?).FontSize = 12;',1:numel(h))

for ip = 1:numel(h1)
  h1(ip).Position(1) = 0.3;
  h1(ip).Position(3) = 0.4;
end
%view([0 0 1])

if 1 % Reset camlight
  %%
  hlight = findobj(h,'type','light');
  delete(hlight)
  for ip = 1:nt
    %camlight(h(ip),45,45);
    camlight(h(ip),45,10);
    camlight(h(ip),45,190);
  end
end

if 0 % Put axis at origin, missing labels, i.e. v_L (km/s)
  %%
  set(gca,'Visible','off')
  %c_eval('h(?).XTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).YTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).ZTick = [-10:2:10]*1e3;',1:numel(h))
  hold(gca,'on')
  % DRAW AXIS LINEs
  plot3(get(gca,'XLim'),[0 0],[0 0],'k');
  plot3([0 0],[0 0],get(gca,'ZLim'),'k');
  plot3([0 0],get(gca,'YLim'),[0 0],'k');
  % GET TICKS
  X=get(gca,'Xtick');
  Y=get(gca,'Ytick');
  Z=get(gca,'Ztick');
  % GET LABELS
  XL=get(gca,'XtickLabel');
  YL=get(gca,'YtickLabel');
  ZL=get(gca,'ZtickLabel');
  % REMOVE TICKS
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  set(gca,'Ztick',[]);
  % GET OFFSETS
  Xoff=diff(get(gca,'XLim'))./30;
  Yoff=diff(get(gca,'YLim'))./30;
  Zoff=diff(get(gca,'ZLim'))./30;
  % DRAW TICKS
  %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
  for i=1:length(X)
     plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
  end;
  for i=1:length(Y)
     plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
  end;
  for i=1:length(Z)
     plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
  end;
  % DRAW LABELS
  text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
  text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
  text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

  xlim_ = get(gca,'XLim');
  ylim_ = get(gca,'YLim');
  zlim_ = get(gca,'ZLim');
  text(gca,xlim_(2),0,0,'L (km/s)')
  text(gca,0,ylim_(2),0,'M (km/s)')
  text(gca,0,0,zlim_(2),'N (km/s)')
  
  hold(gca,'off')
end


fontsize = 16;
legs = {'(a)','(b)','(c)','(d)','(e)'};
ileg = 0;
for ip = 1:numel(h1)
  ileg = ileg + 1;
  irf_legend(h1(ip),legs{ileg},[0.02 0.98],'fontsize',fontsize)
end
for ip = 1:numel(h)
  ileg = ileg + 1;
  irf_legend(h(ip),legs{ileg},[0.02 0.98],'fontsize',fontsize)
end

%% Plot isosurfaces of distribution, PDist, cleaned
units = irf_units;
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');

%times = EpochTT(['2017-07-11T22:34:01.300000000Z';
%  '2017-07-11T22:34:01.800000000Z';...
%  '2017-07-11T22:34:02.300000000Z';...
%  '2017-07-11T22:34:02.800000000Z';...
%  '2017-07-11T22:34:03.300000000Z']);
%times = times + 0.090;
times = times(1:3);
nt = times.length;

vint_red = [-1000 1000];
vint = [3500 5000];
vint = [3000 5000] + 500;
eint = [50 90];
eint = [50 90];
vlim_red = [-10 10];
eint = units.me*(vint*1e3).^2/2/units.eV;
c_eval('PDist = ePDist?;',ic)
c_eval('dmpaB = dmpaB?;',ic)
c_eval('dslE = dslE?;',ic)
c_eval('gseE = gseE?;',ic)
c_eval('scPot = scPot?;',ic)


%if exist('h','var'); delete(h)

nrows = 1;
ncols = nt;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;
%nt = 1;

set(gcf,'defaultLegendAutoUpdate','off');

for it = 1:nt % 3D isosurface
    
  hca = h(isub); isub = isub + 1;
  time = times(it);
  elimit = 40;
  elimit = 35;
  elimit = 30;
  pdist = PDist.elim([elimit Inf]).tlim(time+0.499*0.03*[-1 1]);  
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  pdist.data(:,1:5,:,:) = 3e-27;


  scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); Bn = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); En = E/norm(E);
  %E = gseE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); En = E/norm(E);
  Eperp = cross(Bn,cross(E,Bn)); Eperpn = Eperp\norm(Eperp);
  ExB = cross(E,B); ExBn = Eperp\norm(Eperp);
  
  nSmooth = 3;
  
  iso_values = [0.5e-27, 6e-27];
  %iso_values = [1.5e-30 22e-30];
  iso_values = iso_values(2); 
  set(hca,'colororder',[1 0.5 0.5; 0.5 1 0.5]) 
  set(hca,'colororder',[0.2 0.5 1; 1 0.5 0.5; 0.5 1 0.5]) 
  %Tgse = [Lgse;Mgse;Ngse];
  c_eval('Tdsl = [tsLdsl?.resample(pdist).data;tsMdsl?.resample(pdist).data;tsNdsl?.resample(pdist).data];',ic)
  hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);
  %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',Tdsl);
  
  hs.Patch(1).FaceAlpha = 0.9;
  %hs.Patch(1).FaceAlpha = 0.3;
  %hs.Patch(2).FaceAlpha = 0.9;
  
  legs = arrayfun(@(x)sprintf('%g',x),hs.Values,'UniformOutput',false);
  hleg = legend(hs.Patch,legs,'location','northeast');
  hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
  
  
  if 1 % plot E, ExB, B    
    hold(hca,'on')
    qmax = 10000; 
    nE = E/norm(E);
    nEperp = Eperp/norm(Eperp);
    nB = B/norm(B);
    nExB = ExB/norm(ExB);
    colors_quiv_B = [0 0 0];
    colors_quiv_E = [1 0 0];
    colors_quiv_Eperp = [1 0.5 0.5];
    colors_quiv_ExB = [0 0.5 1]; 

    if 1
      nE = nE*hs.Rotation';
      nExB = nExB*hs.Rotation';
      nB = nB*hs.Rotation';
      nEperp = nEperp*hs.Rotation';
    end
    
    fontsize = 13;
    fontweight = 'bold';
    quiver3(hca,qmax*nB(1)*-1,qmax*nB(2)*-1,qmax*nB(3)*-1,qmax*nB(1)*2,qmax*nB(2)*2,qmax*nB(3)*2,'color',colors_quiv_B,'linewidth',2)
    %text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,'B','color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,sprintf('B=%.0fnT',norm(B)),'color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*nE(1)*-1,qmax*nE(2)*-1,qmax*nE(3)*-1,qmax*nE(1)*2,qmax*nE(2)*2,qmax*nE(3)*2,'color',colors_quiv_E,'linewidth',2)
    %text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,'E','color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,sprintf('E=%.0fmV/m',norm(E)),'color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nEperp(1)*-1,qmax*nEperp(2)*-1,qmax*nEperp(3)*-1,qmax*nEperp(1)*2,qmax*nEperp(2)*2,qmax*nEperp(3)*2,'color',colors_quiv_Eperp,'linewidth',2)
    %text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,'E_\perp','color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,['E_\perp' sprintf('=%.0fmV/m',norm(Eperp))],'color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nExB(1)*-1,qmax*nExB(2)*-1,qmax*nExB(3)*-1,qmax*nExB(1)*2,qmax*nExB(2)*2,qmax*nExB(3)*2,'color',colors_quiv_ExB,'linewidth',2)    
    text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,'E\times{B}','color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)
    %text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,sprintf('E\times{B}=%.0fmV/m',norm(Eperp)),'color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')
    %set(hca,'colororder',colors_quiv)
    %irf_legend(hca,{'B','E','E_\perp','E\times{B}'},[0.98 0.98])
  end
  if 0 % plot LMN (obs in GSE!!)
    %%
     if 1
      L_ = L*hs.Rotation';
      M_ = M*hs.Rotation';
      N_ = N*hs.Rotation';
      
     end

    hold(hca,'on')
    colors_quiv_LMN = [0.5 0.5 0.5];
    quiver3(hca,qmax*L_(1)*-1,qmax*L_(2)*-1,qmax*L_(3)*-1,qmax*L_(1)*2,qmax*L_(2)*2,qmax*L_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*L_(1)*1,qmax*L_(2)*1,qmax*L_(3)*1,'L','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*M_(1)*-1,qmax*M_(2)*-1,qmax*M_(3)*-1,qmax*M_(1)*2,qmax*M_(2)*2,qmax*M_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*M_(1)*1,qmax*M_(2)*1,qmax*M_(3)*1,'M','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*N_(1)*-1,qmax*N_(2)*-1,qmax*N_(3)*-1,qmax*N_(1)*2,qmax*N_(2)*2,qmax*N_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*N_(1)*1,qmax*N_(2)*1,qmax*N_(3)*1,'N','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')

  end
end

for ip = 1:numel(h)
  h(ip).XLabel.String = 'v_L (km/s)';
  h(ip).YLabel.String = 'v_M (km/s)';
  h(ip).ZLabel.String = 'v_N (km/s)';
  vlim = [-9 9]*1e3;
  vlim = 10*[-1 1]*1e3;
  h(ip).XLim = vlim;
  h(ip).YLim = vlim;
  h(ip).ZLim = vlim;
  h(ip).XTick = [-10:5:10]*1e3;
  h(ip).YTick = [-10:5:10]*1e3;
  h(ip).ZTick = [-10:5:10]*1e3;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).ZGrid = 'on';
  h(ip).Box = 'on';
end

hlinks = linkprop(h,{'view'});
c_eval('h(?).FontSize = 12;',1:numel(h))

%view([0 0 1])

if 1 % Reset camlight
  %%
  hlight = findobj(h,'type','light');
  delete(hlight)
  for ip = 1:nt
    %camlight(h(ip),45,45);
    camlight(h(ip),45,10);
    camlight(h(ip),45,190);
  end
end

if 0 % Put axis at origin, missing labels, i.e. v_L (km/s)
  %%
  set(gca,'Visible','off')
  %c_eval('h(?).XTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).YTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).ZTick = [-10:2:10]*1e3;',1:numel(h))
  hold(gca,'on')
  % DRAW AXIS LINEs
  plot3(get(gca,'XLim'),[0 0],[0 0],'k');
  plot3([0 0],[0 0],get(gca,'ZLim'),'k');
  plot3([0 0],get(gca,'YLim'),[0 0],'k');
  % GET TICKS
  X=get(gca,'Xtick');
  Y=get(gca,'Ytick');
  Z=get(gca,'Ztick');
  % GET LABELS
  XL=get(gca,'XtickLabel');
  YL=get(gca,'YtickLabel');
  ZL=get(gca,'ZtickLabel');
  % REMOVE TICKS
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  set(gca,'Ztick',[]);
  % GET OFFSETS
  Xoff=diff(get(gca,'XLim'))./30;
  Yoff=diff(get(gca,'YLim'))./30;
  Zoff=diff(get(gca,'ZLim'))./30;
  % DRAW TICKS
  %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
  for i=1:length(X)
     plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
  end;
  for i=1:length(Y)
     plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
  end;
  for i=1:length(Z)
     plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
  end;
  % DRAW LABELS
  text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
  text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
  text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

  xlim_ = get(gca,'XLim');
  ylim_ = get(gca,'YLim');
  zlim_ = get(gca,'ZLim');
  text(gca,xlim_(2),0,0,'L (km/s)')
  text(gca,0,ylim_(2),0,'M (km/s)')
  text(gca,0,0,zlim_(2),'N (km/s)')
  
  hold(gca,'off')
end
%% Plot isosurfaces of distribution, PDist, cleaned, also with TSeries of locations, 3 sp at the same time
units = irf_units;
ic = [2 3];
ic_str = strrep(num2str(ic),' ','');
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z';...
             '2015-10-16T13:07:02.310Z'];
times = irf_time(times_utc,'utc>EpochTT');

%times = EpochTT(['2017-07-11T22:34:01.300000000Z';
%  '2017-07-11T22:34:01.800000000Z';...
%  '2017-07-11T22:34:02.300000000Z';...
%  '2017-07-11T22:34:02.800000000Z';...
%  '2017-07-11T22:34:03.300000000Z']);
times = times + 0*0.030;
times = times(1:5);

tint_plot = [times.start times.stop] + [-0.5 0.5];
nt = times.length;


colors = pic_colors('matlab');
colors_dists = [0.99 0.6 0.5; 0.6 0.8 0.6];
colors = colors([2 5],:);


%if exist('h','var'); delete(h)

nrows = numel(ic);
ncols = nt;
%h = setup_subplots(nrows,ncols,'vertical');
%[h1,h] = initialize_combined_plot('topbottom',2,nrows,ncols,0.5,'horizontal');
[h1,h] = initialize_combined_plot('topbottom',2,nrows,ncols,0.3,'vertical');

isub = 1;

hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors(ic_str))
plB = cell([]);
c_eval('lmnB? = gseB?*Tgse''; plB{end+1} = lmnB?.x;',ic)
irf_plot(hca,plB,'comp')
hca.YLabel.String = 'B_{LMN} (nT)';  hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',colors)
irf_legend(hca,{'L','M','N'},[0.98 0.98])

if 0
hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',colors)
plE = cell([]);
c_eval('lmnE? = gseE?*Tgse''; plE{end+1} = lmnE?.abs;',ic)
irf_plot(hca,plE,'comp')
set(hca,'ColorOrder',colors)
hold(hca,'on')
%irf_plot(hca,plE.resample(PDist,'mean'))
hold(hca,'off')
hca.YLabel.String = 'E_{LMN} (mV/m)';  hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',colors)
%irf_legend(hca,{'L','M','N'},[0.98 0.98])'
end

if 1
hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',colors)
plE = cell([]);
c_eval('lmnE? = gseVe?*Tgse''; plE{end+1} = lmnE?.abs;',ic)
irf_plot(hca,plE,'comp')
set(hca,'ColorOrder',colors)
hold(hca,'on')
%irf_plot(hca,plE.resample(PDist,'mean'))
hold(hca,'off')
hca.YLabel.String = 'v_{e} (mV/m)';  hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',colors)
irf_legend(hca,{'L','M','N'},[0.98 0.98])'
end

irf_zoom(h1,'x',tint_plot)
irf_zoom(h1,'y')

colors = get(gca,'colororder');
for it = 1:nt % mark locations of distributions
  for ip = 1:numel(h1)
    c_eval('tcenter = ePDist?.tlim(times(it)+0.03*0.49*[-1 1]).time;',ic)
    time_int = tcenter+0.03*0.5*[-1 1];
    hmark = irf_pl_mark(h1(ip),time_int',colors(it,:),'facealpha',0.5);
  end
end


isub = 1;
%nt = 1;

set(gcf,'defaultLegendAutoUpdate','off');

for it = 1:nt % 3D isosurface    
  ic_count = 0;
  for iic = ic
    ic_count = ic_count + 1;
    color = mms_colors(string(iic));
    color = colors_dists(ic_count,:);

    vint_red = [-1000 1000];
    vint = [3500 5000];
    vint = [3000 5000] + 500;
    eint = [50 90];
    eint = [50 90];
    vlim_red = [-10 10];
    eint = units.me*(vint*1e3).^2/2/units.eV;
    c_eval('PDist = ePDist?;',iic)
    c_eval('dmpaB = dmpaB?;',iic)
    c_eval('dslE = dslE?;',iic)
    c_eval('gseE = gseE?;',iic)
    c_eval('scPot = scPot?;',iic)

    hca = h(isub); isub = isub + 1;
    time = times(it);
    time = time;
    elimit = 40;
    elimit = 35;
    elimit = 00;
    %pdist = PDist.elim([elimit Inf]).tlim(time+0.499*0.03*[-1 1]);  
    pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
    pdist.data(:,1:5,:,:) = 3e-27;
  
    tint_tmp = time+0.5*0.03*[-1 1];
    scP = scPot.tlim(tint_tmp); scP = mean(scP.data,1);
    B = dmpaB.tlim(tint_tmp); B = mean(B.data,1); Bn = B/norm(B);
    E = dslE.tlim(tint_tmp); E = mean(E.data,1); En = E/norm(E);
    %E = gseE.tlim(tint_tmp); E = mean(E.data,1); En = E/norm(E);
    Eperp = cross(Bn,cross(E,Bn)); Eperpn = Eperp\norm(Eperp);
    ExB = cross(E,B); ExBn = Eperp\norm(Eperp);
    
    nSmooth = 3;
    
    iso_values = [0.5e-27, 6e-27];
    iso_values = [0.5e-27, 6e-27];
    %iso_values = [1.5e-30 22e-30];
    iso_values = iso_values(2); 
    set(hca,'colororder',[1 0.5 0.5; 0.5 1 0.5]) 
    set(hca,'colororder',colors) 
    %set(hca,'colororder',colors(it,:)) 
    %Tgse = [Lgse;Mgse;Ngse];
    c_eval('Tdsl = [tsLdsl?.resample(pdist).data;tsMdsl?.resample(pdist).data;tsNdsl?.resample(pdist).data];',ic)
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',Tdsl);
    
    hs.Patch(1).FaceColor = color.^1;
    hs.Patch(1).FaceAlpha = 1.0;  
    %hs.Patch(1).FaceAlpha = 0.3;
    %hs.Patch(2).FaceAlpha = 0.9;
    
    legs = arrayfun(@(x)sprintf('%g',x),hs.Values,'UniformOutput',false);
    hleg = legend(hs.Patch,legs,'location','northeast');
    hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
    
    
    if 1 % plot E, ExB, B    
      hold(hca,'on')
      qmax = 10000; 
      nE = E/norm(E);
      nEperp = Eperp/norm(Eperp);
      nB = B/norm(B);
      nExB = ExB/norm(ExB);
      colors_quiv_B = [0 0 0];
      colors_quiv_E = [1 0 0];
      colors_quiv_Eperp = [1 0.5 0.5];
      colors_quiv_ExB = [0 0.5 1]; 
  
      if 1
        nE = nE*hs.Rotation';
        nExB = nExB*hs.Rotation';
        nB = nB*hs.Rotation';
        nEperp = nEperp*hs.Rotation';
      end
      
      fontsize = 13;
      fontweight = 'bold';
      quiver3(hca,qmax*nB(1)*-1,qmax*nB(2)*-1,qmax*nB(3)*-1,qmax*nB(1)*2,qmax*nB(2)*2,qmax*nB(3)*2,'color',colors_quiv_B,'linewidth',2)
      %text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,'B','color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)
      text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,sprintf('B=%.0fnT',norm(B)),'color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)
  
      quiver3(hca,qmax*nE(1)*-1,qmax*nE(2)*-1,qmax*nE(3)*-1,qmax*nE(1)*2,qmax*nE(2)*2,qmax*nE(3)*2,'color',colors_quiv_E,'linewidth',2)
      %text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,'E','color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
      text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,sprintf('E=%.0fmV/m',norm(E)),'color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
      
      quiver3(hca,qmax*nEperp(1)*-1,qmax*nEperp(2)*-1,qmax*nEperp(3)*-1,qmax*nEperp(1)*2,qmax*nEperp(2)*2,qmax*nEperp(3)*2,'color',colors_quiv_Eperp,'linewidth',2)
      %text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,'E_\perp','color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
      text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,['E_\perp' sprintf('=%.0fmV/m',norm(Eperp))],'color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
      
      quiver3(hca,qmax*nExB(1)*-1,qmax*nExB(2)*-1,qmax*nExB(3)*-1,qmax*nExB(1)*2,qmax*nExB(2)*2,qmax*nExB(3)*2,'color',colors_quiv_ExB,'linewidth',2)    
      text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,'E\times{B}','color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)
      %text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,sprintf('E\times{B}=%.0fmV/m',norm(Eperp)),'color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)
  
      hold(hca,'off')
      %set(hca,'colororder',colors_quiv)
      %irf_legend(hca,{'B','E','E_\perp','E\times{B}'},[0.98 0.98])
    end
    if 0 % plot LMN (obs in GSE!!)
    %%
     if 1
      L_ = L*hs.Rotation';
      M_ = M*hs.Rotation';
      N_ = N*hs.Rotation';
      
     end

    hold(hca,'on')
    colors_quiv_LMN = [0.5 0.5 0.5];
    quiver3(hca,qmax*L_(1)*-1,qmax*L_(2)*-1,qmax*L_(3)*-1,qmax*L_(1)*2,qmax*L_(2)*2,qmax*L_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*L_(1)*1,qmax*L_(2)*1,qmax*L_(3)*1,'L','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*M_(1)*-1,qmax*M_(2)*-1,qmax*M_(3)*-1,qmax*M_(1)*2,qmax*M_(2)*2,qmax*M_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*M_(1)*1,qmax*M_(2)*1,qmax*M_(3)*1,'M','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*N_(1)*-1,qmax*N_(2)*-1,qmax*N_(3)*-1,qmax*N_(1)*2,qmax*N_(2)*2,qmax*N_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*N_(1)*1,qmax*N_(2)*1,qmax*N_(3)*1,'N','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')

  end

  end
end
%%
for ip = 1:numel(h)
  h(ip).XLabel.String = 'v_L (km/s)';
  h(ip).YLabel.String = 'v_M (km/s)';
  h(ip).ZLabel.String = 'v_N (km/s)';
  vlim = [-8 8]*1e3;
  vlim = 7*[-1 1]*1e3;
  h(ip).XLim = vlim;
  h(ip).YLim = vlim;
  h(ip).ZLim = vlim;
  h(ip).XTick = [-10:5:10]*1e3;
  h(ip).YTick = [-10:5:10]*1e3;
  h(ip).ZTick = [-10:5:10]*1e3;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).ZGrid = 'on';
  h(ip).Box = 'on';
end

hlinks = linkprop(h,{'view'});
c_eval('h(?).FontSize = 12;',1:numel(h))

%view([0 0 1])
%%
if 1 % Reset camlight
  %%
  hlight = findobj(h,'type','light');
  delete(hlight)
  for ip = 1:numel(h)
    %camlight(h(ip),45,45);
    camlight(h(ip),45,10);
    camlight(h(ip),45,190);
  end
end

if 0 % Put axis at origin, missing labels, i.e. v_L (km/s)
  %%
  set(gca,'Visible','off')
  %c_eval('h(?).XTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).YTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).ZTick = [-10:2:10]*1e3;',1:numel(h))
  hold(gca,'on')
  % DRAW AXIS LINEs
  plot3(get(gca,'XLim'),[0 0],[0 0],'k');
  plot3([0 0],[0 0],get(gca,'ZLim'),'k');
  plot3([0 0],get(gca,'YLim'),[0 0],'k');
  % GET TICKS
  X=get(gca,'Xtick');
  Y=get(gca,'Ytick');
  Z=get(gca,'Ztick');
  % GET LABELS
  XL=get(gca,'XtickLabel');
  YL=get(gca,'YtickLabel');
  ZL=get(gca,'ZtickLabel');
  % REMOVE TICKS
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  set(gca,'Ztick',[]);
  % GET OFFSETS
  Xoff=diff(get(gca,'XLim'))./30;
  Yoff=diff(get(gca,'YLim'))./30;
  Zoff=diff(get(gca,'ZLim'))./30;
  % DRAW TICKS
  %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
  for i=1:length(X)
     plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
  end;
  for i=1:length(Y)
     plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
  end;
  for i=1:length(Z)
     plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
  end;
  % DRAW LABELS
  text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
  text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
  text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

  xlim_ = get(gca,'XLim');
  ylim_ = get(gca,'YLim');
  zlim_ = get(gca,'ZLim');
  text(gca,xlim_(2),0,0,'L (km/s)')
  text(gca,0,ylim_(2),0,'M (km/s)')
  text(gca,0,0,zlim_(2),'N (km/s)')
  
  hold(gca,'off')
end

