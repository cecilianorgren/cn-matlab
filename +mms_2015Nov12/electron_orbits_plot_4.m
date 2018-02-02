%% Electron orbits of several particles, not slice distributions but pitchangle plot
t_center = irf_time('2015-11-12T07:19:21.175000012Z','utc>EpochTT'); 
tintObs = t_center + 1*[-1 1];
%tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
%tintObs = tintObs;
CS_normal_velocity = 70; % km/s

ic = 1;

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsCurvB = mvaCurvB.resample(obsB).tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE.resample(obsB);'...
'obsE = mvaE?_new.tlim(tintObs); obsE = obsE.resample(obsB);'... %'obsE = mvaE?_new.tlim(tintObs) - mvaEht.resample(mvaE?_new).tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;

% Colors
colors = mms_colors('xyz');

% Save old saveParticle (prevent from overwriting) and call this one something else
% call the one produced in this script saveParticle_4
conserve_saveParticle = 1;
if conserve_saveParticle && exist('saveParticle','var') && numel(saveParticle) > 1000
  saveParticle_save = saveParticle;
end

% Model parameters
mms_2015Nov12.Bmodel;

T = 0.2; % integration time


% Plot, 4 panels, incl B
colors = mms_colors('1234');
% Set up plot
nRows = 6;
nCols = 1;
units = irf_units;

clear h;
if 0
  h(1) = subplot(nRows,nCols,1); 
  h(2) = subplot(nRows,nCols,2); 
  h(3) = subplot(nRows,nCols,3); 
  h(4) = subplot(nRows,nCols,4); 
  h(5) = subplot(nRows,nCols,[5 6]); 
  
  h(1).Position(2) = h(1).Position(2)-h(1).Position(4);
  h(1).Position(4) = h(1).Position(4)*1.5;

  dy = 0.000;
  h(2).Position(2) = h(1).Position(2)-h(2).Position(4);
  h(3).Position(2) = h(2).Position(2)-h(3).Position(4);
  h(4).Position(2) = h(3).Position(2)-h(4).Position(4);
  h(5).Position(2) = h(4).Position(2)-h(5).Position(4)-0.08;
elseif 1
  h(1) = subplot(nRows,nCols,1); 
  h(2) = subplot(nRows,nCols,2); 
  h(3) = subplot(nRows,nCols,3); 
  h(4) = subplot(nRows,nCols,[4 5 6]); 
  
  h(1).Position(4) = h(1).Position(4)*1.2;  
  h(2).Position(4) = h(2).Position(4)*1.2;
  h(3).Position(4) = h(3).Position(4)*1.2;
  h(4).Position(4) = h(4).Position(4)*1.2;
  
  h(1).Position(2) = h(1).Position(2)+0.1;
  h(2).Position(2) = h(2).Position(2)+0.1;
  h(3).Position(2) = h(3).Position(2)+0.1; 
  
  h(1).Position(2) = h(1).Position(2)-h(1).Position(4);  
  h(2).Position(2) = h(1).Position(2)-h(2).Position(4);
  h(3).Position(2) = h(2).Position(2)-h(3).Position(4);  
  h(4).Position(2) = h(3).Position(2)-h(4).Position(4);  
else
  h(1) = subplot(nRows,nCols,1); 
  h(2) = subplot(nRows,nCols,2); 
  h(3) = subplot(nRows,nCols,[3 4 5 6]); 
  
  h(1).Position(4) = h(1).Position(4)*1.5;  
  h(2).Position(4) = h(2).Position(4)*1.5;
  h(3).Position(4) = h(3).Position(4)*0.8;
  
  h(1).Position(2) = h(1).Position(2)+0.08;
  h(3).Position(4) = h(3).Position(4)+0.08; 
  
  h(1).Position(2) = h(1).Position(2)-h(1).Position(4);  
  h(2).Position(2) = h(1).Position(2)-h(2).Position(4);
  h(3).Position(2) = h(2).Position(2)-h(3).Position(4);  
end

isub = 1;
if 1 % magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Curvature of magnetic field, obs+mod
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsCurvB.data,'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  modfac = 0.1;
  lineMod = plot(hca,zObs,modCurvB(:,1)*modfac,'--','color',B_colors(1,:));
  plot(hca,zObs,modCurvB(:,2)*modfac,'--','color',B_colors(2,:))
  plot(hca,zObs,modCurvB(:,3)*modfac,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'curb B (1/km)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30]; 
end
if 1 % Electric field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 0 % Time: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end
particle_set = 5;
mms_2015Nov12.electron_orbits_many;
colors = mms_colors('1234');

if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 400];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end
for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  z = saveParticle{iP}.r(:,3); % N
  pa = saveParticle{iP}.pa; % pitch angles    
  colors = mms_colors('1234');
  step = 2;
  hold(hca,'on');
  plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
  hold(hca,'off') 
end
if 0 % density of many particles, starting at +N
  %nP = 200;
  hca = h(isub); isub = isub + 1;
  particle_set = 62;
  mms_2015Nov12.electron_orbits_many;
  %nt(nt==0) = NaN;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt'))
  %colorbar
  shading(hca,'flat')
  view(hca,[0 0 1])
  set(hca,'clim',[0 2000])  
  hca.YLim = [0 180];
  hca.YTick = [45 90 135];
end
if 0 % density of many particles, starting at -N
  %nP = 200;
  hca = h(isub); isub = isub + 1;
  particle_set = 61;
  mms_2015Nov12.electron_orbits_many;
  %nt(nt==0) = NaN;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt'))
  %colorbar
  shading(hca,'flat')
  view(hca,[0 0 1])
  set(hca,'clim',[0 2000])  
  hca.YLim = [0 180];
  hca.YTick = [45 90 135];
end
particle_set = 5;
mms_2015Nov12.electron_orbits_many;
colors = mms_colors('1234');

isub = numel(h);
maxx = 0;
for iP = 1:numel(saveParticle) % Electron trajectories
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);
  hca = h(isub);
  
  if 0
    hl(iP) = plot(hca,x*1e-3,z*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,x(1)*1e-3,z(1)*1e-3,'go',...
              x(end)*1e-3,z(end)*1e-3,'rx'); % plot in km's
  ylabel(hca,'N (km)'); hca.YLim = [-30 30];
  xlabel(hca,'L (km)')  
  else
    hl(iP) =  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
             z(end)*1e-3,x(end)*1e-3,'rx') % plot in km's
  ylabel(hca,'L (km)')
  xlabel(hca,'N (km)'); hca.XLim = [-30 30];
  hca.XLim = h(2).XLim;
  hca.YLim = [30 ceil(maxx*1e-3/10)*10];
  end
  if iP == numel(saveParticle); hold(hca,'off'); end  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
legend(hl,labels{:},'location','northwest')

if conserve_saveParticle
  saveParticle_4 = saveParticle;
  saveParticle = saveParticle_save;
end

% add magnetic field
hold(hca,'on')
yellow = mms_colors('matlab'); yellow = yellow(3,:);
qscale = 0.5;
  plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3,'color',yellow);
  hold(hca,'on')
  nq = 30;
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)    
  dL = 20; % km
  plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3+dL,'color',yellow); 
  quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3 + dL,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)     

c_eval('h(?).XTick = [];',1:3)
c_eval('h(?).XLim = [-30 30];',1:4)
c_eval('h(?).XGrid = ''on'';',1:4)
c_eval('h(?).FontSize = 16;',1:4)
c_eval('h(?).Position(3) = h(2).Position(3);',1:4);
%h(3).Position(2) = h(2).Position(2)-h(2).Position(4)-dy;
%h(4).Position(2) = h(3).Position(2)-h(3).Position(4)-dy;
%h(5).Position(2) = h(4).Position(2)-h(4).Position(4)-dy;
irf_plot_axis_align
c_eval('h(?).Position(1) = 0.15; h(?).Position(3) = 0.7;',1:4)


%% 3D plot of particles only, with magnetic field
figure(76);
hca = subplot(1,1,1);
qscale = 0.5; % for magnetic field quivers
nq = 30;
toplot = fix(linspace(1,size(xyz,1),nq));
Bcolor = mms_colors('matlab'); Bcolor = Bcolor(3,:);

firstPl = 1;
for iP = 1:numel(saveParticle) % Electron trajectories
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);  
  
  switch iP
    case {1 2}
      dm = 0;
    case {3 4}
      dm = 100;
  end
  hl(iP) =  plot3(hca,z*1e-3,-y*1e-3 + dm,x*1e-3,'color',colors(iP,:)); 
  
  if firstPl == 1; firstPl = 0; hold(hca,'on'); axis(hca,'equal'); end
  
  % plot magnetic field lines and quivers 
  if any(iP == [3 4])
  dM = dm+00; % km
  dL = 0;
  dN = 0;
  plot3(hca,xyz(:,3)*1e-3 - dN,-xyz(:,2)*1e-3 + dM,xyz(:,1)*1e-3 - dL,'color',Bcolor);      
  quiver3(hca,...
     xyz(toplot,3)*1e-3 - dN,...
    -xyz(toplot,2)*1e-3 + dM,...
     xyz(toplot,1)*1e-3 - dL,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)
  end
  if any(iP == [1 2])
  dM = 120; % km
  dL = -20;
  dN = 0;
  plot3(hca,xyz(:,3)*1e-3 - dN,-xyz(:,2)*1e-3 + dM,xyz(:,1)*1e-3 - dL,'color',Bcolor);      
  quiver3(hca,...
     xyz(toplot,3)*1e-3 - dN,...
    -xyz(toplot,2)*1e-3 + dM,...
     xyz(toplot,1)*1e-3 - dL,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)
   
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
  
  plot3(hca,z(1)*1e-3,-y(1)*1e-3,x(1)*1e-3,'go',...
            z(end)*1e-3,-y(end)*1e-3,x(end)*1e-3,'rx') % plot in km's
  end
  %hca.XLim = [-30 30];
  %hca.XLim = h(2).XLim;
  %hca.YLim = [30 ceil(maxx*1e-3/10)*10];  
  if iP == numel(saveParticle); hold(hca,'off'); end  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.ZGrid = 'on';
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
legend(hl,labels{:},'location','northwest')

%%
h(2).XTick = [];

h(1).XGrid = 'on';
%h(1).XAxisLocation= 'top';
c_eval('h(?).XLim = h(3).XLim;',[1 4 5]);
h(4).YGrid = 'off';
h(5).YGrid = 'off';
h(5).XLabel.String = 'N (km)';
h(1).YLim = [-13 13];
h(1).Title.String = '';

%c_eval('h(?).CLim = [0 1000];',4:5)
%c_eval('colormap(h(?),''jet'');',4:5)
c_eval('h(?).Box = ''on'';',4:5)


h(end).XLim(1) = 0;
h(end).YGrid = 'on';
h(end).Position(2) = 0.22;

hylab = ylabel(h(2),{'Pitchangle (\circ)'},'interpreter','tex');
hylab.Position(2) = -0.2;

fig = gcf;
hcb = findobj(fig.Children,'type','colorbar');
hcb.Position(1) = h(2).Position(1) + h(2).Position(3)+0.01;
hcb.Position(2) = h(3).Position(2);
hcb.Position(4) = h(2).Position(4)*2;
hcb.Position(3) = 0.025;
hcb.FontSize = 10;
hcb.YLabel.FontSize = 10;

hcb = colorbar('peer',h(4));
% hcb.Position(1) = h(4).Position(1) + h(4).Position(3)+0.01;
% hcb.Position(2) = h(5).Position(2);
% hcb.Position(4) = h(4).Position(4)*2;
 hcb.Position(3) = 0.025;
% hcb.FontSize = 10;
% hcb.YLabel.FontSize = 10;
hcb.YLabel.String = 'Counts';

%cmap = []

hcb = colorbar('peer',h(5));
% hcb.Position(1) = h(4).Position(1) + h(4).Position(3)+0.01;
% hcb.Position(2) = h(5).Position(2);
% hcb.Position(4) = h(4).Position(4)*2;
 hcb.Position(3) = 0.025;
% hcb.FontSize = 10;
% hcb.YLabel.FontSize = 10;
hcb.YLabel.String = 'Counts';


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
legshift = 0; % the two sc configuration plots

for ii = 1:6
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii).FontSize = 12;  
  h(ii).YLabel.FontSize = 12;  
end

set(gcf, 'InvertHardCopy', 'off');
set(gcf,'paperpositionmode','auto');
set(gcf,'color','white');

%% Electron orbits of several particles, not slice distributions but pitchangle plot, 3 planes for particle trajectories
close
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
tintObs = tintObs;
CS_normal_velocity = 70; % km/s

ic = 1;

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsCurvB = mvaCurvB.resample(obsB).tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE.resample(obsB);'...
'obsE = mvaE?_new.tlim(tintObs); obsE = obsE.resample(obsB);'... %'obsE = mvaE?_new.tlim(tintObs) - mvaEht.resample(mvaE?_new).tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;

% Colors
colors = mms_colors('xyz');

% Model parameters
mms_2015Nov12.Bmodel;

T = 0.2; % integration time


% Plot, 4 panels, incl B
scrsz = get(groot,'ScreenSize');
figurePosition = scrsz;
figurePosition(3) = 0.7*scrsz(3);
figure('position',figurePosition)

colors = mms_colors('1234');
% Set up plot
nRows = 6;
nCols = 2;
units = irf_units;

clear h;

  % left column
  h(1) = subplot(nRows,nCols,1); 
  h(2) = subplot(nRows,nCols,3); 
  h(3) = subplot(nRows,nCols,5); 
  h(4) = subplot(nRows,nCols,[7 9 11]); 
  % right column
  h(5) = subplot(nRows,nCols,[2 4 6]); 
  h(6) = subplot(nRows,nCols,[8 10 12]); 
  
  if 1
    h(1).Position(4) = h(1).Position(4)*1.2;  
    h(2).Position(4) = h(2).Position(4)*1.2;
    h(3).Position(4) = h(3).Position(4)*1.2;
    h(4).Position(4) = h(4).Position(4)*1.2;
    h(6).Position(4) = h(6).Position(4)*1.2;

    h(1).Position(2) = h(1).Position(2)+0.1;
    h(2).Position(2) = h(2).Position(2)+0.1;
    h(3).Position(2) = h(3).Position(2)+0.1;     

    h(1).Position(2) = h(1).Position(2)-h(1).Position(4);  
    h(2).Position(2) = h(1).Position(2)-h(2).Position(4);
    h(3).Position(2) = h(2).Position(2)-h(3).Position(4);  
    h(4).Position(2) = h(3).Position(2)-h(4).Position(4);  
    
    h(6).Position(2) = h(4).Position(2);
    h(5).Position(2) = h(3).Position(2);
    h(5).Position(4) = h(3).Position(4)*3;
    
    h(5).Position(3) = 0.4;h(5).Position(3)*1.2;
    h(6).Position(3) = 0.4;h(6).Position(3)*1.2;
    
  end

isub = 1;
if 1 % magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Curvature of magnetic field, obs+mod
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsCurvB.data,'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  modfac = 0.1;
  lineMod = plot(hca,zObs,modCurvB(:,1)*modfac,'--','color',B_colors(1,:));
  plot(hca,zObs,modCurvB(:,2)*modfac,'--','color',B_colors(2,:))
  plot(hca,zObs,modCurvB(:,3)*modfac,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'curb B (1/km)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30]; 
end
if 1 % Electric field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 0 % Time: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  %pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end

particle_set = 5;
mms_2015Nov12.electron_orbits_many;
colors = mms_colors('1234');

for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  z = saveParticle{iP}.r(:,3); % N
  pa = saveParticle{iP}.pa; % pitch angles    
  colors = mms_colors('1234');
  step = 2;
  hold(hca,'on');
  plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
  hold(hca,'off') 
end




maxx = 0;
for iP = 1:numel(saveParticle) % Electron trajectories
  isub = 4;
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);
  hca = h(isub);
  
  if 1
    for is = isub + [0 1 2];
      hca = h(is);
      hl(iP) =  plot3(hca,x*1e-3,y*1e-3,z*1e-3,'color',colors(iP,:));
      if is == 5; hleg(iP) = hl(iP); end
      if iP == 1; hold(hca,'on'); end
      plot3(hca,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go',...
                x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx') % plot start and end in km's
      xlabel(hca,'L (km)')
      ylabel(hca,'M (km)')
      zlabel(hca,'N (km)');
      %hca.ZLim = [-30 30];
      %hca.YLim = [30 ceil(maxx*1e-3/10)*10];

      if iP == numel(saveParticle); hold(hca,'off'); end  
    end
  else 
    % NL
    hca = h(isub); isub = isub + 1;
    hl(iP) =  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
             z(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'N (km)'); hca.XLim = [-30 30];
    hca.XLim = h(2).XLim;
    hca.YLim = [30 ceil(maxx*1e-3/10)*10];

    if iP == numel(saveParticle); hold(hca,'off'); end  

    % MN
    hca = h(isub); isub = isub + 1;
    hl(iP) =  plot(hca,-y*1e-3,-z*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,-y(1)*1e-3,-z(1)*1e-3,'go',...
             -y(end)*1e-3,-z(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'N (km)'); %hca.XLim = [-30 30];
    %hca.XLim = h(2).XLim;
    %hca.YLim = [30 ceil(maxx*1e-3/10)*10];

    if iP == numel(saveParticle); hold(hca,'off'); end  

    % ML
    hca = h(isub); isub = isub + 1;
               plot(hca,-y*1e-3,x*1e-3,'color',colors(iP,:));
    %hl(iP) =  plot(hca,-y*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,-y(1)*1e-3,x(1)*1e-3,'go',...
             -y(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'-M (km)');
    %hca.XLim = h(2).XLim;
    %hca.YLim = [30 ceil(maxx*1e-3/10)*10];    
  end
  
  
  if iP == numel(saveParticle); hold(hca,'off'); end    
  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
legend(hleg,labels{:},'location','northwest')

% add magnetic field
for iP = 4:6
  if 1 % 3    
    nq = 20;
    qscale  = 0.5;
    
    toplot = fix(linspace(1,size(xyz,1),nq));
    toplot = [1 fix(size(xyz,1)*0.85)];
    hca = h(iP);
    hold(hca,'on')
    dL = 0;
    dM = 75; % km
    dN = 0;
    plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor);      
    quiver3(hca,...
       xyz(toplot,1)*1e-3 - dN,...
       xyz(toplot,2)*1e-3 - dM,...
       xyz(toplot,3)*1e-3 - dL,...
       Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       qscale,'color',Bcolor)
        
    dL = 0;
    dM = 10; % km
    dN = 0;
    plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor);      
    quiver3(hca,...
       xyz(toplot,1)*1e-3 - dN,...
       xyz(toplot,2)*1e-3 - dM,...
       xyz(toplot,3)*1e-3 - dL,...
       Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       qscale,'color',Bcolor)
     hold(hca,'off')
  else % 2D
    hca = h(iP);
    hold(hca,'on')
    yellow = mms_colors('matlab'); yellow = yellow(3,:);
    qscale = 0.5;
    plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3,'color',yellow);
    hold(hca,'on')
    nq = 30;
    toplot = fix(linspace(1,size(xyz,1),nq));
    quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3,...
       Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)    
    dL = 20; % km
    plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3+dL,'color',yellow); 
    quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3 + dL,...
       Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)     
  end
end

c_eval('h(?).XTick = [-20 0 20];',1:2)
c_eval('h(?).XTick = [];',3)
c_eval('h(?).XLim = [-30 30];',1:3)
%c_eval('h(?).XLim = [0 115];',4:6)
%c_eval('h(?).YLim = [-50 110];',4:6)
%c_eval('h(?).ZLim = [-32 32];',4:6)
c_eval('h(?).XTick = -200:20:200;',4:6)
c_eval('h(?).YTick = -200:20:200;',4:6)
c_eval('h(?).ZTick = -200:20:200;',4:6)

irf_plot_axis_align(h(1:3))

c_eval('h(?).XGrid = ''on'';',1:6)
c_eval('h(?).YGrid = ''on'';',4:6)
c_eval('h(?).ZGrid = ''on'';',4:6)
c_eval('h(?).FontSize = 14;',1:6)
c_eval('h(?).Box = ''on'';',4:6)
c_eval('h(?).Position(3) = h(2).Position(3);',1:6);
h(5).Position(3) = 0.4;
h(6).Position(3) = 0.4;

h(5).YTickLabel = [];
%h(3).Position(2) = h(2).Position(2)-h(2).Position(4)-dy;
%h(4).Position(2) = h(3).Position(2)-h(3).Position(4)-dy;
%h(5).Position(2) = h(4).Position(2)-h(4).Position(4)-dy;

view(h(4),[0 -1 0]); h(4).ZDir = 'reverse'; camroll(h(4),90); h(4).XAxisLocation = 'top';
view(h(5),[1 0 0]); 
view(h(6),[0 0 -1]); camroll(h(6),90); h(6).XAxisLocation = 'top';



