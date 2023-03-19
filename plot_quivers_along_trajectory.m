% From 
L = [0.85 0.52 0.10];
M = [-0.49 0.80 0.01];
N = [-0.08 -0.06 0.94];

L = [1 0 0];
M = [0 1 0];
N = [0 0 1];
if 0
L = [0.9578   -0.2873         0];
M = [0.2873    0.9578         0];
N = [     0         0    1.0000];
end
rot = [L;M;N];

ic = 1;
colors = pic_colors('matlab');

tint_plot = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');
tint_plot = tint_plot + [-1 1];
dt = 0.2; % s
timeline = tint_plot(1):dt:tint_plot(end);
  
% Define progression of line to plot from
c_eval('vQ = gseVExB?*rot'';',ic)
vQ = vQ.resample(timeline);
xQ = -cumsum(vQ.x.data)*dt; % km/s*s = km
yQ = -cumsum(vQ.y.data)*dt; % km/s*s = km
zQ = -cumsum(vQ.z.data)*dt; % km/s*s = km

% Define quantities to plot as quivers
Q = cell(0,0);
Qleg = cell(0,0);
if 1
c_eval('Q(end+1) = {gseB?*rot''};  Qleg(end+1) = {''B''};',ic)
c_eval('Q(end+1) = {gseE?*rot''};  Qleg(end+1) = {''E''};',ic)
c_eval('Q(end+1) = {gseVe?*rot''};  Qleg(end+1) = {''Ve''};',ic)
c_eval('Q(end+1) = {gseVi?*rot''};  Qleg(end+1) = {''Vi''};',ic)
c_eval('Q(end+1) = {gseJ?*rot''};  Qleg(end+1) = {''J''};',ic)

else % different matlab versions?
c_eval('Q{end+1} = gseB?*rot'';  leg{end+1} = ''B'';',ic)
c_eval('Q{end+1} = gseE?*rot'';  leg{end+1} = ''E'';',ic)
c_eval('Q{end+1} = gseVe?*rot''; leg{end+1} = ''Ve'';',ic)
c_eval('Q{end+1} = gseVi?*rot''; leg{end+1} = ''Vi'';',ic)
c_eval('Q{end+1} = gseJ?*rot'';  leg{end+1} = ''J'';',ic)
end
Q_ = cell(0,0);
c_eval('Q_{?} = Q{?}.resample(timeline);',1:numel(Q))

% Define quantities to plot as symbols, e.g. circles or squares
S = cell(0,0);
Sleg = cell(0,0);
Symb = cell(0,0);
c_eval('stmp = gseE?perp*rot''; S{end+1} = stmp.z; Sleg{end+1} = ''E''; Symb{end+1} = ''o'';',ic)

hca = subplot(1,1,1);
axis(hca,'equal')
hca.NextPlot = 'add'; % add or replace
%set(hca,'ColorOrder',colors)
for iQ = 1:numel(Q)  
  pQ = Q{iQ}.resample(timeline).data;  
  quiver3(hca,xQ,yQ,zQ,pQ(:,1),pQ(:,2),pQ(:,3),'color',colors(iQ,:),'linewidth',1) 
end
for iS = 1:numel(S)
  pS = S{iS}.resample(timeline).data;  
  scatter3(hca,xQ,yQ,zQ,abs(pS(:,1))*10,[1 0 0],'filled') 
end
hleg = legend(hca,Qleg,'box','off','location','best');

hold(hca,'off')
%hca.NextPlot = 'replace'; % add or replace

hca.XLabel.String = 'x';
hca.YLabel.String = 'y';
hca.ZLabel.String = 'z';
hca.XDir = 'reverse';
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.ZGrid = 'on';
hca.Box = 'on';
hca.FontSize = 16;
