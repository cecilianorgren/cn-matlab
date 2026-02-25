%% Figure 4: For paper, 1 sc, maybe 2 sc, horizontal reversed time positioning
scList = 1;
ic = scList(1);
tintZoom = irf.tint('2015-10-16T10:33:28.00Z/2015-10-16T10:33:31.00Z');

% Initialize plot
nRows = 3;
nCols = 4;

clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end


isub = 1;

% Plot ve time series
hca = h(isub); isub = isub + 1;
c_eval('irf_plot(hca,{mvaVe?par,mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)
irf_legend(hca,{'v_{||}','v_{\perp,L}','v_{\perp,M}','v_{\perp,N}'},[0.02 0.1])
hca.YLabel.String = {'v_e','(km/s)'};
%irf_legend(hca,{irf_ssub('MMS ?',ic)},[0.1 0.95])
hca.Title.String = irf_ssub('MMS ?',ic);

hca = h(isub); isub = isub + 1;
c_eval('irf_plot(hca,{fce?*2*pi},''comp'');',ic)
hca.YLabel.String = {'\omega_{ce}','(1/s)'};
hca.Title.String = irf_ssub('MMS ?',ic);

hca = h(isub); isub = isub + 1;
c_eval('irf_plot(hca,{Le?,re?},''comp'');',ic)
irf_legend(hca,{'L_{e}','\rho_{e}'},[0.02 0.1])
hca.YLabel.String = {'Length','(km)'};
hca.Title.String = irf_ssub('MMS ?',ic);

nTS = isub - 1;

irf_zoom(h(1:nTS),'x',tintZoom)
irf_zoom(h(1:nTS),'y')
irf_plot_axis_align

% Decide times for projection plots
indTimes = 1207:1215; zdists = [1 1 2 3 4 5 6 7 8];
indTimes = 1206:1214; zdists = [1 1 1 2 3 4 5 6 7]+0;

for iicc = 1:numel(scList)
  ic = scList(iicc);
  c_eval('dist = ePDist?;',ic)
  times = dist.time;
  indTime = indTimes;%indTimes{ic};
  % Projection coordinate system
  csys = 5;
  switch csys
    case 1 % z: B, x: N, y: zxx    
      c_eval('z = gseB?.resample(ePDist?);',ic); 
      z = -z/z.abs;
      tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
      tsL = irf.ts_vec_xyz(z.time,repmat(L,z.length,1));

      x = cross(z,cross(tsN,z));    
      y1 = cross(z,x);
      y2 = cross(z,cross(tsL,z));    
      y = y2;
      vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
    case 2
      x = hatB0;
      y = hatExB0;
      z = cross(x,y);
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
    case 5 % z: B, x: N, y: zxx    
      tref = irf_time('2015-10-16T10:33:20.266Z','utc>epochTT');
      c_eval('z = gseB?.resample(ePDist?);',ic); 
      z = -z/z.abs;
      %tsN = irf.ts_vec_xyz(z.time,repmat(N,z.length,1));
      c_eval('tsL = gseVe?.resample(tref);',ic); tsL = tsL/tsL.abs;
      tsL = irf.ts_vec_xyz(z.time,repmat(tsL.data,z.length,1));
      y = cross(z,cross(-tsL,z))
      x = cross(z,y);    
      %y1 = cross(z,x);
      %y2 = cross(z,cross(tsL,z));    
      %y = y2;
      vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
  end
  X = x;
  Y = y;
  Z = z;

  c_eval('dist = ePDist?;',ic)
  c_eval('scpot = scPot?;',ic)
  vlim = 15*1e3;
  elevlim = 10;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  haveYLabel = 0;
  for ii = 1:numel(indTime) % plot mms1, plane: NxB, N 
    
    x = X(indTime(ii)).data;
    y = Y(indTime(ii)).data;
    z = Z(indTime(ii)).data;

    time = times(indTime(ii));
    timeUTC = time.utc;  
    hmark = irf_pl_mark(h(1:nTS),[time.epochUnix]),mms_colors(irf_ssub('?',ic));

    % Get mean vectors
    c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
    hatVe0 = double(irf_norm(Ve0));    
    c_eval('B0 = gseB?.resample(time).data;',ic); 
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = gseE?.resample(time).data;',ic); 
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);
   % vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};


    hca = h(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot);    
    %hca.Title.String = timeUTC(15:23);
    hca.Title.String = '';
    ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    ht.VerticalAlignment = 'top';
    ht.HorizontalAlignment = 'left';
    ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';
    
    % Mark velocities and different distributions on plot
    hold(hca,'on')
    % Core 
    angles = 0:360;
    v0 = 5.2; % 10^3 km/s
    x0 = v0*cosd(angles);
    y0 = v0*sind(angles);
    plot(hca,x0,y0,'w');
    
     % Outside accelerated edge
    angles = 0:360;
    v1 = 8.5; % 10^3 km/s
    x1 = v1*cosd(angles);
    y1 = v1*sind(angles);
    plot(hca,x1,y1,'k');
    
    % Get z from difference in circles above
    % E = kz; 
    % (ek/m)z^2 = (e/m)Ez = (e/m)phi
    ekmz2 = (v1*1e6)^2-(v0*1e6)^2; % (m/s)^2
    phi = ekmz2*units.me/units.e;
    
    
    c_eval('oce = fce?.resample(time).data*2*pi;',ic)
    c_eval('ExB = gseVExB?.resample(time).data;',ic)
    c_eval('E = gseE?perp.resample(time).abs;',ic)
    c_eval('B = gseB?.resample(time).abs;',ic)
    
    
    % electron penetration depth lim
    v0 = 100:100:10000; % km/s
    ze = 4/oce*(v0+norm(ExB));
    
    
    % 'vperp' lim
    zdist = 6; % km
    zdist = zdists(ii);
    vz = -10e3:1e3:10e3;
    vy = vz.^2/oce/zdist-norm(ExB)-0.25*oce*zdist;
    plot(hca,vz*1e-3,vy*1e-3,'k')
    
    % Add some info
    ht = text(hca.XLim(1),hca.YLim(2),['z = ' num2str(zdist,'%g') ' km'],'color',[1 1 1]);
    ht.VerticalAlignment = 'top';
    ht.HorizontalAlignment = 'right';
    ht.FontSize = 12;    
    
    ht = text(hca.XLim(1),hca.YLim(1),['v_e = ' num2str(norm(Ve0),'%.0f') ' km/s'],'color',[1 1 1]);
    ht.VerticalAlignment = 'bottom';
    ht.HorizontalAlignment = 'right';
    ht.FontSize = 12;    
    
    hold(hca,'off')
    
    % Calculate moments
    
    
    hca.XLabel.String = 'N_{\perp}';
    hca.YLabel.String = 'N_{\perp}\times B';
       
    if 0
      hca = h2(isub); isub = isub + 1; 
      %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      mms.plot_projection(hca,dist,'tint',time,'xyz',[-y;-z;-x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);
      titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
      hca.Title.String = titleStr;
      hca.Title.String = '';
      colormap(hca,strCMap)
      hca.XDir = 'reverse';

      hca.XLabel.String = 'N_{\perp}';
      hca.YLabel.String = 'B';
    end
  end   
end


hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);

cmap = irf_colormap('space');
cmap = 'jet';
for ii = 1:nRows*nCols
  h(ii).FontSize = 12;
  colormap(h(11),cmap)
end
