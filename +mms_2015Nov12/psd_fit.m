ic = 3;

n = [5]*1e6;
m = [0];
t = [50]*1e-3;
vd = [0.3];
d = [2];
a1 = [1];
a2 = [1];
toPlot = [1];


% Real distribution
time = irf_time('2015-11-12T07:19:21.162Z','utc>epochtt');
c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)
tInd = ePitch.time.tlim(time+0.015*[-1 1]);

nRows = 1; nCols = 2;
for  ii = 1:nRows*nCols
  h(ii) = subplot(nRows,nCols,ii);
end

isub = 1;

hca = h(isub); isub = isub + 1;
colors = hca.ColorOrder;
if 1
toPlot = [1];
whamp.plot_f(hca,n(toPlot),m(toPlot),t(toPlot),vd(toPlot),d(toPlot),a1(toPlot),a2(toPlot),'pitchangles',[0 90 180],'PSDvsE','km/s');
hca.YScale = 'log';
hca.XScale = 'log';
hca.XLim = [1e1 2e3];
hca.YLim = [1e-2  1e5];
hca.XTick = [1e-1 1e0 1e1 1e2 1e3];
end


hold(hca,'on')
%hl = plot(hca,ePitch(905).depend{1},[ePitch(905).data(:,:,1); mean(ePitch(905).data(:,:,9:10),3); ePitch(905).data(:,:,18)],'*');
hl = plot(hca,ePitch(905).depend{1},[ePitch(905).data(:,:,1)],'*');
hl.Color = colors(1,:);
hl = plot(hca,ePitch(905).depend{1},[mean(ePitch(905).data(:,:,9:10),3)],'*');
hl.Color = colors(2,:);
hl = plot(hca,ePitch(905).depend{1},[ePitch(905).data(:,:,18)],'*');
hl.Color = colors(3,:);
hca.XScale = 'log';
hca.YScale = 'log';
hold(hca,'off')


hca = h(isub); isub = isub + 1;
plot(hca,ePitch(905).depend{1},squeeze(ePitch(905).data(:,:,:)))
hca.XScale = 'log';
hca.YScale = 'log';
hca.XLim = [1e1 2e3];
hca.YLim = [1e-2  1e5];
hca.XTick = [1e-1 1e0 1e1 1e2 1e3];

 
