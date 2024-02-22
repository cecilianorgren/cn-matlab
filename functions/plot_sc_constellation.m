function varargout = plot_sc_constellation(varargin)

% Check for axes
[ax,args,nargs] = irf.axescheck(varargin{:}); 
if isempty(ax), ax = gca; end
hca = ax;

R1 = args{1};
R2 = args{2};
R3 = args{3};
R4 = args{4};
      
% Spacecraft constellation
R0 = (R1 + R2.resample(R1) + R3.resample(R1) + R4.resample(R1))/4;
c_eval('r? = R?-R0;',1:4)
c_eval('r? = r?.resample(R0).data;',1:4)
%c_eval('dr(?,!) = dr? - dr!;',1:4,1:4)
r = [r1;r2;r3;r4];
rmax = max(r);
rmin = min(r);


% Plotting
symbols = {'o','s','^','<'};
colors(1,:) = [0.2 0.2 0.2];

colors = [0.2000    0.2000    0.2000;
          0.9000    0.2000         0;
               0    0.8000         0;
          0.1000    0.4000    0.9000];
linewidth = 2;
markersize = 20;

holdon = 0;
for isc = 1:4
  plot3(hca,r(isc,1),r(isc,2),r(isc,3),'s','markersize',markersize,'linewidth',linewidth,'color',colors(isc,:));
  if not(holdon)
    hold(hca,'on')
  end
end

hca.XLim = [-15 15];
hca.YLim = [-15 15];
hca.ZLim = [-15 15];

hca.XLim = [rmin(1) rmax(1)]+[-1 1]*2;
hca.YLim = [rmin(2) rmax(2)]+[-1 1]*2;
hca.ZLim = [rmin(3) rmax(3)]+[-1 1]*2;

xmin = hca.XLim(2);
ymin = hca.YLim(2);
zmin = hca.ZLim(1);
for isc = 1:4
  plot3(hca,xmin*[1 1],r(isc,2),r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),ymin*[1 1],r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),r(isc,2),zmin*[1 1],'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
end

for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],'linewidth',1,'color',[0 0 0]);
  end
end
for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[xmin xmin],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[ymin ymin],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[zmin zmin],':','linewidth',1,'color',[1 1 1]*0.7);
  end
end

hold(hca,'off')

hca.XLabel.String = 'x (km)';
hca.YLabel.String = 'y (km)';
hca.ZLabel.String = 'z (km)';

%hca(?).FontSize = 14;',1:numel(h))

hca.Box = 'on';
%axis(hca,'square');
%daspect(hca,[1 1 1])
hca.DataAspectRatio = [1 1 1];

