%% Load data
tint = irf.tint('2015-11-01T15:08:00.00Z/2015-11-01T15:08:15.00Z');
ic = 4;

c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('tic; ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?); toc',ic)
%c_eval('[ePDist?,ePDistError?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

units = irf_units;
c_eval('dslVExB? = cross(dslE?.resample(dmpaB?.time),dmpaB?)/dmpaB?.abs2*1e3;',ic) % km/s
c_eval('oce? = units.e*dmpaB?.abs*1e-9/units.me; oce?.name = ''omega ce''; oce?.units = ''rad/s'';',ic); % Hz

c_eval('dslE?slow = dslE?.resample(ePDist?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(ePDist?);',ic)

%% Figure: Crescent distributions and Bessho fit, 1 single time
% Decide times
time = irf_time('2015-11-01T15:08:06.824Z','utc>epochtt');
time=time+-0*0.03*10;%*3;
d = 4; % km

% Initialize particle distribution plot
nRows = 2;
nCols = 4;
isub = 0;
% Make the ordering vertical, so one can plot different 
% planes one time in one column
doVerticalSorting = 0;
if doVerticalSorting % first fill 1st column, then 2nd column etc.
  for jj = 1:nCols
    for ii = 1:nRows  
      isub = isub + 1;         
      h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
    end
  end
else % like normal subplot(-,-,-): % first fill 1st row, then 2nd row etc.
  for ii = 1:nRows
    for jj = 1:nCols
      isub = isub + 1;         
      h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
    end
  end
end

isub = 1;

c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('VExB = dslVExB?.resample(time).abs.data;',ic)
c_eval('oce = oce?.resample(time).data;',ic)
%VExB = 0;

% Projection coordinate system
c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
c_eval('hatExB = cross(hatE,hatB);',ic)

c_eval('[out,l,v] = irf_minvar(dmpaB?,''<Bn>=0'');',ic)
L = v(1,:); M = v(2,:); N = v(3,:);
Lperp = cross(hatB,cross(L,hatB)); Lperp = Lperp/norm(Lperp);
Nperp = cross(hatB,cross(N,hatB)); Nperp = Nperp/norm(Nperp);
Mperp = cross(hatB,cross(M,hatB)); Mperp = Mperp/norm(Mperp);

csys = 1;
switch csys  
  case 1 % ExB, 3rd, B    
    y = hatE;
    z = hatB;
    x = cross(z,x)/norm(cross(z,x));
    vlabels = {'E\times B','E','B'};
  case 2
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 3
    
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 4.5];  
  

timeUTC = time.utc;      

% Perpendicular plane    
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';
hca.Title.String = timeUTC(1:23);

% B plane 1
hca = h2(isub); isub = isub + 1; 
xyz = [x;z;-y]; vlabels = {'v_{ExB}','v_B','-v_{E}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% B plane 2
hca = h2(isub); isub = isub + 1;
xyz = [y;z;x]; vlabels = {'v_E','v_B','v_{ExB}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% Perp plane, but with markings
vectors = {hatExB,'ExB'; hatE,'E'};
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
hca.Title.String = '';

% Mark velocities and different distributions on plot
hold(hca,'on')

% 'vperp' lim: Bessho2016
zdist = d;
vz = -10e3:1e3:10e3;
%vy = vz.^2/oce/zdist-VExB-0.25*oce*zdist;
vy = @(vz,d) vz.^2/oce/d-VExB-0.25*oce*d;
ang =  -90;-30; % rotate parabola as is necessary
vz_rot = @(vz,vy,ang) vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
vy_rot = @(vz,vy,ang) -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
vz_rot = @(vz,d,ang) vz*1e-3*cosd(ang)+vy(vz,d)*1e-3*sind(ang);
vy_rot = @(vz,d,ang) -vz*1e-3*sind(ang)+vy(vz,d)*1e-3*cosd(ang);

%newvz = vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
%newvy = -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
%h_vcrescent = plot(hca,newvz,newvy,'k');
delete(h_vcrescent)
if 1
  ds = [3 4 5 6];
  id = 0;
  for ii = 1:numel(ds)
    id = id +1;
    h_vcrescent1 = plot(hca,vz_rot(vz,ds(id),ang),vy_rot(vz,3,ang),'color',[0 0 0]); 
  end
  hca.Title.String = {['d = [' num2str(ds,'%4.1f') '] km'],'furthest in -> furthest out'};
elseif 1
  h_vcrescent1 = plot(hca,vz_rot(vz,3,ang),vy_rot(vz,3,ang),'color',[     0    0.4470    0.7410]);
  h_vcrescent2 = plot(hca,vz_rot(vz,4,ang),vy_rot(vz,4,ang),'color',[0.4940    0.1840    0.5560]);%[0.8500    0.3250    0.0980]);
  h_vcrescent3 = plot(hca,vz_rot(vz,5,ang),vy_rot(vz,5,ang),'color',[0.9290    0.6940    0.1250]);[0.4660    0.6740    0.1880]
  h_vcrescent4 = plot(hca,vz_rot(vz,6,ang),vy_rot(vz,6,ang),'color',[0.8500    0.3250    0.0980]);
  hleg=legend([h_vcrescent1 h_vcrescent2 h_vcrescent3 h_vcrescent4],{'d = 3','d = 4','d = 5','d = 6'},'location','eastoutside');
  hleg.Position = hleg.Position + [0.07 0 0 0];
elseif 1
  h_vcrescent1 = plot(hca,vz_rot(vz,vy(vz,3),ang),vy_rot(vz,vy(vz,3),ang),'color',[     0    0.4470    0.7410]);
  h_vcrescent2 = plot(hca,vz_rot(vz,vy(vz,4),ang),vy_rot(vz,vy(vz,4),ang),'color',[0.4940    0.1840    0.5560]);%[0.8500    0.3250    0.0980]);
  h_vcrescent3 = plot(hca,vz_rot(vz,vy(vz,5),ang),vy_rot(vz,vy(vz,5),ang),'color',[0.9290    0.6940    0.1250]);[0.4660    0.6740    0.1880];
  h_vcrescent4 = plot(hca,vz_rot(vz,vy(vz,6),ang),vy_rot(vz,vy(vz,6),ang),'color',[0.8500    0.3250    0.0980]);
  hleg=legend([h_vcrescent1 h_vcrescent2 h_vcrescent3 h_vcrescent4],{'d = 3','d = 4','d = 5','d = 6'},'location','eastoutside');
  hleg.Position = hleg.Position + [0.07 0 0 0];
else
  h_vcrescent = plot(hca,vz_rot(vz,vy(vz,3),ang),vy_rot(vz,vy(vz,3),ang),...
                     hca,vz_rot(vz,vy(vz,4),ang),vy_rot(vz,vy(vz,4),ang),...
                     hca,vz_rot(vz,vy(vz,5),ang),vy_rot(vz,vy(vz,5),ang));
  hleg=legend([h_vcrescent],{'d = 3','d = 4','d = 5'},'location','eastoutside');
  hleg.Position = hleg.Position + [0.07 0 0 0];    
end

% Add some info
%hca.Title.String = ['d = ' num2str(zdist,'%g') ' km'];
hold(hca,'off')

% Perpendicular plane, fixed coordinate system
hca = h2(isub); isub = isub + 1; 
xyz = [L;M;N]; vlabels = {'v_{L}','v_M','v_N'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';
hca.Title.String = timeUTC(1:23);

% B plane 1
hca = h2(isub); isub = isub + 1; 
xyz = [L;N;-M]; vlabels = {'v_{L}','v_N','-v_{M}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% B plane 2
hca = h2(isub); isub = isub + 1;
xyz = [M;N;L]; vlabels = {'v_M','v_N','v_L'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% Perp plane, but with markings
hca = h2(isub); isub = isub + 1; 
% Mark velocities and coordinate systems
if 0
  hold(hca,'on')
  vectors = [L;M;N;Lperp;Mperp;Nperp;hatB;hatE;hatExB];
  quiver3(hca,vectors(:,1)*0,vectors(:,1)*0,vectors(:,1)*0,...
              vectors(:,1),vectors(:,2),vectors(:,3))
else
  quiver3(hca,0,0,0,L(1),L(2),L(3),'k'); hold(hca,'on')
  quiver3(hca,0,0,0,M(1),M(2),M(3),'k')
  quiver3(hca,0,0,0,N(1),N(2),N(3),'k')
  
  quiver3(hca,0,0,0,Lperp(1),Lperp(2),Lperp(3),'r')
  quiver3(hca,0,0,0,Mperp(1),Mperp(2),Mperp(3),'r')
  quiver3(hca,0,0,0,Nperp(1),Nperp(2),Nperp(3),'r')
  
  quiver3(hca,0,0,0,hatB(1),hatB(2),hatB(3),'g')
  quiver3(hca,0,0,0,hatE(1),hatE(2),hatE(3),'y')
  quiver3(hca,0,0,0,hatExB(1),hatExB(2),hatExB(3),'c')
end
hold(hca,'off')
axis(hca,'equal')
hca.XGrid = 'on'; hca.YGrid = 'on'; hca.ZGrid = 'on';

% Add some info
%hca.Title.String = ['d = ' num2str(zdist,'%g') ' km'];

 

colormap(strCMap)
% print
% cn.print(['crescent_fit_' irf_time(time,'epochtt>utc_yyyymmddTHHMMSS_mmm') '_d' num2str(d,'%g') 'km'],'time',time)


%% Figure 1: Crescent distributions and
% run first mms_2015Oct16.partial_moments
ic = 1;
zdists = {sort([7 6.5 5.5 4.5 3.5 3 1.5 0.1 0.1])};
indTimes = {sort([1213 1212 1211 1210 1209 1208 1207 1206 1205])-5*3*0};

step = 2;
zdists = {zdists{ic}(1:step:end)}; indTimes = {indTimes{ic}(1:step:end)};

% Initialize particle distribution plot
nRows = 4;
nCols = 5;
isub = 0;
% initialize plot, make the ordering vertical, so one can plot different 
% planes one time in one column
for jj = 1:nCols
  for ii = 1:nRows  
    isub = isub + 1;         
    h2(isub) = subplot(nRows,nCols,jj+(ii-1)*nCols);    
  end
end

% Decide times
%indTimes = {[1214 1212 1210 1208 1206-10]+1,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
%indTimes = {[1214 1212 1210 1208 1206 1204],[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};

isub = 1;

c_eval('scpot = scPot?;',ic)
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
times = dist.time;
indTime = indTimes{ic};

% Projection coordinate system
csys = 1;
switch csys
  case 1 % z: B, x: N, y: zxx    
    c_eval('z = gseB?.resample(ePDist?);',ic); 
    z = -z/z.abs;
    tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
    tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

    x = cross(z,cross(tsN,z));    
    y1 = cross(z,x);
    y2 = cross(z,cross(tsL,z));    
    y = y1;
    vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  case 2 % ExB, 3rd, B
    x = hatExB0;
    y = cross(z,x);
    z = hatB0;
    vlabels = {'B','E\times B','B\times(E\times B)'};
  case 3
    x = [1 0 0];
    y = [0 1 0];
    z = [0 0 1];
    vlabels = {'X','Y','Z'};
  case 4
    x = L;
    y = M;
    z = N;
    vlabels = {'L','M','N'};
end
X = x;
Y = y;
Z = z;

vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

haveYLabel = 0;
for ii = 1:nCols % loop over times
  % plot mms1, plane: NxB, N 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;

  time = times(indTime(ii));
  timeUTC = time.utc;      

  % Get mean vectors
  c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
  hatVe0 = double(irf_norm(Ve0));    
  c_eval('B0 = gseB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = gseE?.resample(time).data;',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
 % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};
   
  if nRows > 0 % Perpendicular plane    
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'v_{N\times B} (10^3 km/s)';                    
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
  end
  if nRows > 1 % B plane
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[x;z;-y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
      hca.Title.String = '';
      %hca.Title.String = timeUTC(15:23);
      %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'left';
      %ht.FontSize = 12;
      colormap(hca,strCMap)
      hca.XDir = 'reverse';
      if ii == 1;          
        hca.YLabel.String = 'v_{B} (10^3 km/s)';          
      else
        hca.YLabel.String = '';
        hca.YTickLabel = ''; 
      end    
      hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
    end
  if nRows > 2 % B plane
    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);        
    hca.Title.String = '';
    %hca.Title.String = timeUTC(15:23);
    %ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    %ht.VerticalAlignment = 'top';
    %ht.HorizontalAlignment = 'left';
    %ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      hca.YLabel.String = 'v_{B} (10^3 km/s)';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    %hca.XLabel.String = 'N_{\perp}\times B';   
    hca.XLabel.String = 'v_{N\times B} (10^3 km/s)';   
  end
  if nRows > 3 % Perp plane, but with markings and data of partial moments
    c_eval('neCore = ne?Core.resample(time).data;',ic)
    c_eval('neCrescent = ne?Crescent.resample(time).data;',ic)
    c_eval('veCore = ve?Core.resample(time).data;',ic)
    c_eval('veCrescent = ve?Crescent.resample(time).data;',ic)
    vectors = {hatExB0,'ExB'; E0,'E'}; % veCrescent,['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'];

    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vectors',vectors);        
    hca.Title.String = '';

    % Mark velocities and different distributions on plot
    hold(hca,'on')

    c_eval('oce = fce?.resample(time).data*2*pi;',ic)
    c_eval('ExB = gseVExB?.resample(time).data;',ic)
    c_eval('E = gseE?perp.resample(time).abs;',ic)
    c_eval('B = gseB?.resample(time).abs;',ic)     

    c_eval('v0 = minV?.resample(time).data;',ic) % km/s
    % Core 
    angles = 0:360;
    %v0 = 5.5*1e3; % 10^3 km/s
    x0 = v0*cosd(angles);
    y0 = v0*sind(angles);
    h_vcirc = plot(hca,x0*1e-3,y0*1e-3,'k');

     % Outside accelerated edge
    angles = 0:360;
    v1 = 8.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');
    
    angles = 0:360;
    v1 = 7.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    %plot(hca,x1,y1,'k');

    % 'vperp' lim
    zdist = 6; % km
    zdist = zdists{ic}(ii);
    vz = -10e3:1e3:10e3;
    vy = vz.^2/oce/zdist-norm(ExB)-0.25*oce*zdist;
    ang =  -15; -90+acosd(x*irf_norm(Ve0)');
    newvz = vz*1e-3*cosd(ang)+vy*1e-3*sind(ang);
    newvy = -vz*1e-3*sind(ang)+vy*1e-3*cosd(ang);
    h_vcrescent = plot(hca,newvz,newvy,'k');
    if ii == 1
      delete(h_vcrescent)
    end

    % Add some info
    if 0 % separate spots
      ht = text(hca.XLim(2),hca.YLim(2),['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],'color',[1 1 1]);
      ht.VerticalAlignment = 'top';
      ht.HorizontalAlignment = 'left';
      ht.FontSize = 12; 

      htd = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      htd.VerticalAlignment = 'top';
      htd.HorizontalAlignment = 'right';
      htd.FontSize = 12;          

      ht = text(hca.XLim(1),hca.YLim(1),['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s'],'color',[1 1 1]);
      ht.VerticalAlignment = 'bottom';
      ht.HorizontalAlignment = 'right';
      ht.FontSize = 12;        
    elseif 1 % all gathered

      if 0
        ht = text(hca.XLim(1),hca.YLim(2),{['n_c = ' num2str(neCrescent,'%.2g') ' cm^{-3}'],['v_e = ' num2str(norm(veCrescent),'%.0f') ' km/s']},'color',[1 1 1]);
        ht.VerticalAlignment = 'top';
        ht.HorizontalAlignment = 'right';
        ht.FontSize = 9; 
      end
      %ht = text(hca.XLim(1),hca.YLim(2),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
      %ht.VerticalAlignment = 'top';
      %ht.HorizontalAlignment = 'right';
      %ht.FontSize = 12;   

      if 0
        htd = text(hca.XLim(1),hca.YLim(1),['d = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
        htd.VerticalAlignment = 'bottom';
        htd.HorizontalAlignment = 'right';
        htd.FontSize = 10;       
        htd.FontWeight = 'bold'; 
      else % put it int he title
        hca.Title.String = ['d = ' num2str(zdist,'%g') ' km'];
        %htd.FontSize = 10;       
        hca.Title.FontWeight = 'light'; 
        if ii ==1
          hca.Title.String = '';
        end
      end
    end

    hold(hca,'off')


    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    if ii == 1;          
      %hca.YLabel.String = 'N_{\perp}\times B';          
      hca.YLabel.String = 'v_{N\times B} (10^3 km/s)';          
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end    
    hca.XLabel.String = 'v_{N_{\perp}} (10^3 km/s)';   
  end
  end   


for ii = 1:(nRows*nCols)
  originalPosition{ii} = h2(ii).Position;  
end
