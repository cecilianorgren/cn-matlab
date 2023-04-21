% Solar rotation

% degrees/day
A = 14.713; % pm 0.0491
B = -2.396; % pm 0.188
C = -1.787; % pm 0.253


w = @(lat) A + B*sind(lat).^2 + C*sind(lat).^4;
lat = 0:90;

fontsize = 14;
colors = pic_colors('matlab');
doMarkKp = 1;
TKp = [9,13.6,27];

nrows = 2;
ncols = 1;
h = gobjects([nrows,ncols]);
h(1) = subplot(nrows,ncols,1);
h(2) = subplot(nrows,ncols,2);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,lat,w(lat))
hca.XLabel.String = 'Solar latitude (degrees)';
irf_legend(hca,{'Rotation frequency','(degrees/day)'}',[0.98 0.98],'fontsize',fontsize,'color',colors(1,:))
%hca.YLabel.String = {'Rotation frequency','(degrees/day)'};
hca.YLim = [10 15];

hca = h(isub); isub = isub + 1;
plot(hca,lat,360./w(lat))
hca.YTick = 24:2:36;
hca.YLim = [24 35];
hca.XLabel.String = 'Solar latitude (degrees)';
%hca.YLabel.String = 'Rotation period (days)';
irf_legend(hca,{'Rotation period','(days)'}',[0.02 0.98],'fontsize',fontsize,'color',colors(1,:))
if doMarkKp
  hold(hca,'on')
  for ii = 1:numel(TKp)
    plot(hca,[0 90],TKp(ii)*[1 1])
  end
  hold(hca,'off')
end

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).XLim = [0 90];',1:numel(h))
c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).YGrid = ''on'';',1:numel(h))

compact_panels(0.03,0.05)