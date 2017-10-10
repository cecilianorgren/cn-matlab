%% Electron orbits of several particles, not slice distributions but pitchangle plot
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
tintObs = tintObs;
CS_normal_velocity = 70; % km/s

ic = 1;

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE.resample(obsB);'...
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
colors = mms_colors('1234');
% Set up plot
nRows = 7;
nCols = 1;
units = irf_units;

h(1) = subplot(nRows,nCols,1); 
h(2) = subplot(nRows,nCols,2); 
h(3) = subplot(nRows,nCols,3); 
h(4) = subplot(nRows,nCols,4); 
h(5) = subplot(nRows,nCols,5); 
h(6) = subplot(nRows,nCols,[6 7]); 

isub = 1;
if 1 % magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,20)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
end
if 1 % Time: ePDist pa low energies
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
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.CLim = h(2).CLim;  
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
%%
for iP = 1:numel(saveParticle) % Electron trajectories
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  
  hca = h(isub);
  
  if 1
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
  end
  if iP == numel(saveParticle); hold(hca,'off'); end  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
legend(hl,labels{:},'location','west')

h(1).Position(2) = h(1).Position(2)-h(1).Position(4);
h(1).Position(4) = h(1).Position(4)*2;



dy = 0.000;
h(2).Position(2) = h(1).Position(2)-h(1).Position(4)*0.5-dy;
h(3).Position(2) = h(2).Position(2)-h(2).Position(4)-dy;
h(4).Position(2) = h(3).Position(2)-h(3).Position(4)-dy;
h(5).Position(2) = h(4).Position(2)-h(4).Position(4)-dy;


h(2).Position(2) = h(1).Position(2)-h(1).Position(4)-dy;
h(3).Position(2) = h(2).Position(2)-h(2).Position(4)-dy;
h(4).Position(2) = h(3).Position(2)-h(3).Position(4)-dy;
h(5).Position(2) = h(4).Position(2)-h(4).Position(4)-dy;

yy = 0.08;
h(2).Position(4) = yy;
h(3).Position(4) = yy;
h(4).Position(4) = yy;
h(5).Position(4) = yy;

h(2).Position(2) = yy;
h(3).Position(2) = yy;
h(4).Position(2) = yy;
h(5).Position(2) = yy;

for ii = 1:6;
  h(ii).Position(3) = 0.7;
end

h(2).XTick = [];

h(1).XGrid = 'on';
%h(1).XAxisLocation= 'top';
c_eval('h(?).XLim = h(3).XLim;',[1 4 5]);
h(4).YGrid = 'off';
h(5).YGrid = 'off';
h(5).XLabel.String = 'N (km)';
h(1).YLim = [-13 13];
h(1).Title.String = '';

c_eval('h(?).CLim = [0 1000];',4:5)
c_eval('colormap(h(?),''jet'');',4:5)
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

cmap = []

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




