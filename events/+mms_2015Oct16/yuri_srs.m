%% Figure 4: For paper, 1 sc, horizontal reversed time positioning, only perpendicular plane
scList = [1];
ic = scList;
nsc = numel(scList);  

% Initialize particle distribution plot
nRows = 1;
nCols = 4;
isub = 1;
clear h2;

for ii = 1:(nRows*nCols);
  h2(isub) = subplot(nRows,nCols,ii); 
  h2(isub).Position(3) = h2(isub).Position(3)*1.2;
  h2(isub).Position(1) = h2(isub).Position(1)-0.03;
  isub = isub + 1; 
end


% Decide times
indTimes = {[1214 1212 1210 1208 1206-10]+4,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
indTimes = {[1214 1212 1210 1208 1206]+0,[1214 1212 1210 1208 1206-2]+1,[1205 1207 1194 1186 1184],[1205 1207 1194 1186 1184]};
indTimes = {[1213 1212 1211 1210]};
isub = 1;
for iicc = 1:numel(scList)
  ic = scList(iicc);
  c_eval('dist = ePDist?;',ic)
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
      y = y2;
      vlabels = {'B\times(N\times B)','B\times(B\times(N\times B))','B'};
    case 2
      c_eval('x = irf_norm(dslE?);',ic)
      z = irf_norm(gseB1);
      y = cross(z,x);
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

  c_eval('dist = ePDist?;',ic)
  vlim = 15*1e3;
  elevlim = 10;
  strCMap = 'jet';
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  haveYLabel = 0;
  for ii = 1:numel(h2) % plot mms1, plane: NxB, N 
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


    hca = h2(isub); isub = isub + 1; 
    %mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim);    
    %hca.Title.String = timeUTC(15:23);
    hca.Title.String = '';
    ht = text(hca.XLim(2),hca.YLim(2),timeUTC(15:23),'color',[1 1 1]);
    ht.VerticalAlignment = 'top';
    ht.HorizontalAlignment = 'left';
    ht.FontSize = 12;
    colormap(hca,strCMap)
    hca.XDir = 'reverse';    


    if ~haveYLabel          
      hca.YLabel.String = 'N_{\perp}\times B';    
      hca.YLabel.String = 'L_{\perp}';    
      haveYLabel = 1;
    else
      hca.YLabel.String = '';
      hca.YTickLabel = ''; 
    end
    %hca.XTickLabel = ''; 
    %hca.XLabel.String = ''; 
    hca.XLabel.String = 'N_{\perp}';    
  end 
  haveYLabel = 0;
  for ii = 1:numel(indTime) % plot mms1, plane: B, ... 
  x = X(indTime(ii)).data;
  y = Y(indTime(ii)).data;
  z = Z(indTime(ii)).data;
    
  time = times(indTime(ii));
  timeUTC = time.utc;    
end 
end

hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCB.Position(1) = hCB.Position(1)+0.08;
%hCB.Position = h2(end).Position;
%hCB.Position(1) = hCB.Position(1)+hCB.Position(3);

%%
for ii = 1:nRows*nCols
  originalPosition{ii} = h2(ii).Position;
end


%hCB.Position = [h2(nCols).Position(1)+xWidth+0.02 h2(4).Position(2) 0.02 h2(4).Position(4)];

cmap = irf_colormap('space');
cmap = 'jet';
for ii = 1:nRows*nCols
  h2(ii).Position(3) = originalPosition{ii}(3)*2.7;
  h2(ii).Position(4) = originalPosition{ii}(4)*1.5;
  h2(ii).Position(2) = originalPosition{ii}(2)-0.05;
  h2(ii).Position(1) = originalPosition{ii}(1)-0.02;
  
  h2(ii).FontSize = 12;
  colormap(h2(11),cmap)
end

for ii = [6:10 16:20]
  h2(ii).Position(2) = h2(ii).Position(2)+0.026;
end
  
for ii = 11:20
  h2(ii).Position(2) = h2(ii).Position(2)-0.05;
end