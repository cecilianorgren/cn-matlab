nA = 360;
nBin = 32;
dA = nA/nBin;
allA = 0:dA:(360-dA);
deg_per_sec = 360/20;

T = 20; % s
tOneBin = 0.150;
dA_from_spinning = 0.150;

colors = pic_colors('matlab');

hca = subplot(1,1,1);
plot(hca,cosd(0:360),sind(0:360),'-k')

hold(hca,'on')
scale = 1/0.9;
iTurn = 0;
for t = 0:tOneBin:(4*tOneBin)  
  iTurn = iTurn + 1;
  color = colors(iTurn,:);

  deltaAngle = t*deg_per_sec;

  scale = scale*0.9;

  plot(hca,scale*cosd(0:360),scale*sind(0:360),'color',color)
  xpatch = scale*[0 sind(allA(1)+deltaAngle) sind(allA(2)+deltaAngle) 0];
  ypatch = scale*[0 cosd(allA(1)+deltaAngle) cosd(allA(2)+deltaAngle) 0];
  patch(hca,xpatch,ypatch,color)
  for iA = 1:nBin
    plot(hca,scale*[0 sind(allA(iA)+deltaAngle)],scale*[0 cosd(allA(iA)+deltaAngle)],'color',color)
  end
  pause(1)
end
hold(hca,'off')

hca.XAxisLocation = "origin";
hca.YAxisLocation = "origin";
hca.XTickLabel = [];
hca.YTickLabel = [];

%%

tt = irf_time('2017-07-11T22:33:50.000Z','utc>EpochTT');
elim = [5000 15000];
nt = 6;
h = setup_subplots(6,1);


for it = 0:(nt-1)
  hca = h(it+1);
  tt_tmp = tt + 0.5*0.150*[-1 1] +it*0.150;
  pdist = iPDist3.elim(elim).tlim(tt_tmp);
  az = pdist.depend{2};   
  pol = pdist.depend{3};
    
  az = [az(1)-11.25/2 az+11.25/2];
  pol = [pol(1)-11.25/2 pol+11.25/2];
  [AZ,POL] = ndgrid(az,pol);

  data = squeeze(mean(pdist.data,2));
  surf(hca,AZ,POL,AZ*0,data)
  view(hca,[0 0 1])
  shading(hca,'flat')
end

drawnow
hlinks = linkprop(h,{'XLim','YLim','CLim','XTick','YTick'});
colormap(pic_colors('candy_gray'))
c_eval('h(?).Layer = "top";',1:numel(h))
c_eval('h(?).Box = "on";',1:numel(h))
c_eval('h(?).FontSize = 12;',1:numel(h))
c_eval('h(?).YLabel.String = {"Polar angle","(deg)"};',1:numel(h))
c_eval('h(?).XLabel.String = "Azimuthal angle (deg)";',numel(h))
h(1).XLim = [0 360]+11.25*[-1 1];
h(1).YLim = [0 180];
h(1).XTick = [0:30:360];
h(1).YTick = [0:30:180];
compact_panels(h,0.01,0.01)
h(1).Title.String = sprintf('%.0f < E < %.0f eV',elim(1),elim(2))

