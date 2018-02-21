cd /Users/cno062/Research/Events/2017-07-11_223323//Whamp

n = [6 14];
m = [0 0];
t = [05 30]*1e-3;
vd = [00 0];
d = [0.001 1];
a1 = [3 0.05];
a2 = [6 1];


h = subplot(2,2,1);
toPlot = [1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar

h = subplot(2,2,2);
toPlot = [2];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar

h = subplot(2,2,3);
toPlot = [2 1];
whamp.plot_f(h,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot));
colorbar


%% run whamp
ttime = irf_time('2017-07-11T22:34:03.000Z','utc>epochtt');

distribution = 1;
switch distribution 
  case 1
    B = 1e-9;
    n = [-0.03 0.03];
    m = [0 1];
    t = [300 8000];
    vd = [1.8 0];
    a = [2 1];  % a1    % anisotropy, Tperp/Tpar 
    d = [1 1];  % d     % loss cone parameter 1, 1 when no loss cone
    b = [0 0];  % b     % loss cone parameter 2, 0 when no loss cone
    toPlot = [1 2];
  case 2
    n = [-24 40];
    m = [0 0];
    t = [23 35];
    vd = [0 0];
    d = [0.2 0.5];
    a1 = [2.25 1.9];
    a2 = [1.1 0];
    toPlot = [1 2];
end


clear Species;
for ii = toPlot
  species.n = n(ii);
  species.m = m(ii);
  species.t = t(ii);
  species.vd = vd(ii);
  species.a = a(ii);
  species.d = d(ii);
  species.b = b(ii);
  Species{ii} = species;
end
  
PlasmaModel.B = 17;
PlasmaModel.Species = Species;

InputParameters.fstart = 0.0001; % Hz,  start frequency
InputParameters.kperp = [0 0.010 1];   % - kp (perpendicular wave vector) as scalar or, [kp_start kp_step kp_end]
InputParameters.kpar = [0.02,0.0050,1.2];
InputParameters.useLog = 0;
Output = whamp.run(PlasmaModel,InputParameters); % run whamp

fce = 0.0279928*1000*PlasmaModel.B; % fce in Hz

% max growth rate, along par direction
iWimax = find(imag(Output.f(:,1))==max(imag(Output.f(:,1))));
colors = mms_colors('matlab');

nRows = 3;
nCols = 2;
for ip = 1:nRows*nCols; h(ip) = subplot(nRows,nCols,ip); end 
isub = 1;

if 1 % Measured electron pitch angle distribution and WHAMP fit
  hca = h(isub); isub = isub +1;  
  whamp.plot_f(hca,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a(toPlot),b(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XLim = [1e1 5e3];
  hca.YLim = [1e-2  1e6];
  hca.XTick = [1e-1 1e0 1e1 1e2 1e3];

  % add real distribution
  time = tint(1); 
  tInd = find(abs(ePitch1.time-time)==min(abs(ePitch1.time-time)));
  timeUTC = time.utc;
  hold(hca,'on')
  unitscale = 1e30; % cm^-6 ->km^-6
  mms_pa = plot(hca,ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,1)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,7)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,13)*unitscale,'+');
  mms_pa(1).Color = hca.ColorOrder(1,:);
  mms_pa(2).Color = hca.ColorOrder(2,:);
  mms_pa(3).Color = hca.ColorOrder(3,:);
  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLim = [1e-32 2e-25]*unitscale;
  hca.XLim = [1e1 3e3];  
  hca.YLabel.String = 'f_e (s^3 km^{-6})';
  hca.XTick = [1e1 1e2 1e3 1e4];
  %hleg = irf_legend(h,{'0';'90';'180'},[0.98 0.98]);
  hca.Title.String = ['Electron distribution'];
  %hold(h,'on')
  %mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])
  %hold(h,'off')
  irf_legend(hca,{{'-  Fit';'+ Measured';['   ' timeUTC(1:10)];['   ' timeUTC(12:23)]}},[0.05 0.4],'color',[0 0 0])
end

isub = isub + 1;
if 1 % Surface plot of dispersion relation
  hca = h(isub); isub = isub +1;  
  
  surf_height = real(Output.f); 
  surf_color = imag(Output.f);
  
  surf_height(surf_height<0) = NaN;
  surf_height(surf_height>5) = NaN;
  
  surf(hca,Output.kperp,Output.kpar,surf_height,surf_color)
  shading(hca,'flat')
  hca.XLabel.String = 'k_{\perp}';
  hca.YLabel.String = 'k_{||}';
  hca.ZLabel.String = 'f_r/f_{ce}';
  hca.XLim = Output.kperp([1 end]);
  hca.YLim = Output.kpar([1 end]);
  hcb = colorbar('peer',hca);
  hcb.YLim = hca.CLim;
  hca.CLim = max(abs(hca.CLim))*[-1 1]*1;
  irf_colormap('poynting')
  hcb.YLabel.String = 'f_i/f_{ce}';
  hca.Title.String = ['Dispersion relation'];
end

if 1 % 2D flat plot of dispersion relation, frequency as contours
  hca = h(isub); isub = isub +1; 
  
  surf_height = real(Output.f); 
  surf_color = imag(Output.f);
  
  surf_height(surf_height<0) = NaN;
  surf_height(surf_height>5) = NaN;
  
  pcolor(hca,Output.kperp,Output.kpar,imag(Output.f))
  hold(hca,'on')
  [c,ha] = contour(hca,Output.kperp,Output.kpar,surf_height,'k');
  clabel(c,ha)
  hold(hca,'on')
  shading(hca,'flat')

  hca.Title.String = ['Dispersion relation'];
  hca.XLabel.String = 'k_{\perp}';
  hca.YLabel.String = 'k_{||}';
  hca.ZLabel.String = 'f/f_{ce}';
  %hca.YScale = 'log';
  %hca.XScale = 'log';
  hca.XLim = Output.kperp([1 end]);
  hca.YLim = Output.kpar([1 end]);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_i/f_{ce}';
  hcb.YLim = hca.CLim;
  hca.CLim = max(abs(hca.CLim))*[-1 1]*0.7;
  cmap = irf_colormap(hca,'poynting');
  %colormap(hca,cmap(1:10:end,:))
  axis(hca,'square')
  %hca.XLim = hca.XLim*0.9;
end

if 1 % 1D plot of parallel dispersion relation
  hca = h(isub); isub = isub +1;
  %set(hca,'ColorOrder',[colors(1:2,:); colors(1:2,:)])
  plot(hca,Output.kpar,(real(Output.f(:,1))),...
           Output.kpar,imag(Output.f(:,1)),...
           Output.kpar(iWimax),real(Output.f(iWimax,1)),'*',...
           Output.kpar(iWimax),imag(Output.f(iWimax,1)),'*')
  irf_legend(hca,{'Real frequency','Growth rate'},[0.05 0.95])
  irf_legend(hca,{sprintf('w_{r,max} = %g',real(Output.f(iWimax,1))),sprintf('w_{i,max} = %g',imag(Output.f(iWimax,1)))},[0.05 0.05])
  hca.Title.String = ['Dispersion relation, k_{\perp}=' num2str(Output.kperp(1),'%g')];
  hca.YLabel.String = 'f/f_{ce}';
  hca.XLabel.String = 'k_{||}';
  hca.YLim = [-0.2 hca.YLim(2)];  
  hca.XLim = [Output.kpar([1 end])];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

if 1 % 1D plot of parallel phase velocity
  hca = h(isub); isub = isub +1;
  plot(hca,Output.kpar,real(Output.f(:,1)./Output.kpar')*2*pi,Output.kpar(iWimax),real(Output.f(iWimax,1)./Output.kpar(iWimax))*2*pi,'*')  
  hca.Title.String = 'Parallel phase velocity';
  hca.YLabel.String = '2\pi f/k_{||}';
  hca.XLabel.String = 'k_{||}';
  hca.YLim = [0 hca.YLim(2)];  
  hca.XLim = [Output.kpar([1 end])];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end


%% never mind about the higher energy ring, jus tfocus on the lower part.
tint = irf.tint('2015-10-16T10:33:44.90Z/2015-10-16T10:33:45.00Z'); % magnetosphere-magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:33:44.90Z/2015-10-16T10:33:44.93Z')+0.03; % magnetosphere-magnetosheath-magnetosphere



n = [15];
m = [0];
t = [26]; % temperature
a = [1.8]; % tperp/tpar
d = [1]; % loss cone parameter 1
b = [0]; % loss cone parameter 2
vd = [0]; % drift speed/thermal speed
toPlot = [1];


    n = [-23 40];
    m = [0 0];
    t = [25 35];
    vd = [0 0];
    d = [0.2 1];%d = 1*[1 1];
    a = [2.2 1.7];
    b = [1.1 0];
    toPlot = [1 2];
    
clear Species species;
for ii = toPlot
  species.n = n(ii);
  species.m = m(ii);
  species.t = t(ii);
  species.a = a(ii);
  species.d = d(ii);
  species.b = b(ii);
  species.vd = vd(ii);
  Species{ii} = species;
end
Species{1};
  
PlasmaModel.B = 17;
PlasmaModel.Species = Species;

InputParameters.fstart = 0.001; % Hz,  start frequency
InputParameters.kperp = [0 0.010 1];   % - kp (perpendicular wave vector) as scalar or, [kp_start kp_step kp_end]
InputParameters.kpar = [0.02,0.0050,1.2];
InputParameters.useLog = 0;
Output = whamp.run(PlasmaModel,InputParameters); % run whamp

fce = 0.0279928*1000*PlasmaModel.B; % fce in Hz

nRows = 2;
nCols = 1;
for ip = 1:nRows*nCols; h(ip) = subplot(nRows,nCols,ip); end 

isub = 1;

hca = h(isub); isub = isub + 1;
[ikmax] = find(imag(Output.f(:,1))==max(imag(Output.f(:,1))));

plot(hca,Output.kpar,[real(Output.f(:,1)) imag(Output.f(:,1))],...
         Output.kpar(ikmax),[real(Output.f(ikmax,1)) imag(Output.f(ikmax,1))],'*')
hca.XLim = [0 1];
hca.YLim = [-0.2 0.4];

if 1 % Measured electron distribution and WHAMP fit
  hca = h(isub); isub = isub +1;
  whamp.plot_f(hca,n(toPlot)*1e6,m(toPlot),t(toPlot)*1e-3,vd(toPlot),d(toPlot),a(toPlot),b(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
  hca.YScale = 'log';
  hca.XScale = 'log';
  hca.XLim = [1e1 5e3];
  hca.YLim = [1e-2  1e6];
  hca.XTick = [1e-1 1e0 1e1 1e2 1e3];

  % add real distribution
  time = tint(1); 
  tInd = find(abs(ePitch1.time-time)==min(abs(ePitch1.time-time)));
  timeUTC = time.utc;
  hold(hca,'on')
  unitscale = 1e30; % cm^-6 ->km^-6
  mms_pa = plot(hca,ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,1)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,7)*unitscale,'+',...
                    ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,13)*unitscale,'+');
  mms_pa(1).Color = hca.ColorOrder(1,:);
  mms_pa(2).Color = hca.ColorOrder(2,:);
  mms_pa(3).Color = hca.ColorOrder(3,:);

  hold(hca,'off')
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLim = [1e-32 2e-25]*unitscale;
  hca.XLim = [1e1 3e3];  
  hca.YLabel.String = 'f_e (s^3 km^{-6})';
  hca.XTick = [1e1 1e2 1e3 1e4];
  %hleg = irf_legend(h,{'0';'90';'180'},[0.98 0.98]);
  hca.Title.String = ['Electron distribution'];
  %hold(h,'on')
  %mms.plot_pitchangles(h,ePDist1,dmpaB1,'tint',tint,'scPot',scPot1,'pitchangle',[0 90 180])
  %hold(h,'off')
  irf_legend(hca,{{'-  Fit';'+ Measured';['   ' timeUTC(1:10)];['   ' timeUTC(12:23)]}},[0.05 0.4],'color',[0 0 0])
end

%%
if 1 % Surface plot of dispersion relation
  hca = h(isub); isub = isub +1;
  surf(hca,Output.kperp,Output.kpar,(real(Output.f)),imag(Output.f))
  shading(hca,'flat')
  hca.XLabel.String = 'k_{\perp}';
  hca.YLabel.String = 'k_{||}';
  hca.ZLabel.String = 'f_r/f_{ce}';
  hca.XLim = Output.kperp([1 end]);
  hca.YLim = Output.kpar([1 end]);
  hcb = colorbar('peer',hca);
  hcb.YLim = hca.CLim;
  hca.CLim = max(abs(hca.CLim))*[-1 1]*1;
  irf_colormap('poynting')
  hcb.YLabel.String = 'f_i/f_{ce}';
  hca.Title.String = ['Dispersion relation'];
end

if 1 % 1D plot of parallel dispersion relation
  hca = h(isub); isub = isub +1;
  plot(hca,Output.kpar,(real(Output.f(:,1))),Output.kpar,imag(Output.f(:,1)))
  irf_legend(hca,{'Real frequency','Growth rate'},[0.05 0.95])
  hca.Title.String = ['Dispersion relation, k_{\perp}=' num2str(Output.kperp(1),'%g')];
  hca.YLabel.String = 'f/f_{ce}';
  hca.XLabel.String = 'k_{||}';
  hca.YLim = [-0.2 hca.YLim(2)];  
  hca.XLim = [Output.kpar([1 end])];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

if 1 % 2D flat plot of parallel dispersion relation, frequency as contours
  hca = h(isub); isub = isub +1; 
  pcolor(hca,Output.kperp,Output.kpar,imag(Output.f))
  hold(hca,'on')
  [c,ha] = contour(hca,Output.kperp,Output.kpar,real(Output.f),'k');
  clabel(c,ha)
  hold(hca,'on')
  shading(hca,'flat')

  hca.Title.String = ['Dispersion relation'];
  hca.XLabel.String = 'k_{\perp}';
  hca.YLabel.String = 'k_{||}';
  hca.ZLabel.String = 'f/f_{ce}';
  %hca.YScale = 'log';
  %hca.XScale = 'log';
  hca.XLim = Output.kperp([1 end]);
  hca.YLim = Output.kpar([1 end]);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_i/f_{ce}';
  hcb.YLim = hca.CLim;
  hca.CLim = max(abs(hca.CLim))*[-1 1]*0.7;
  cmap = irf_colormap(hca,'poynting');
  %colormap(hca,cmap(1:10:end,:))
  axis(hca,'square')
  %hca.XLim = hca.XLim*0.9;
end
