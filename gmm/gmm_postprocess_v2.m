% Save path
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v1/'];
listFiles = dir([pathSave '*iFF*']);
nFiles = numel(listFiles);


out = gmmCrossEventStatistics(listFiles(1:50),'version',{'minbick','chordbick'});
%out = gmmCrossEventStatistics(listFiles(1:10));

%% Figure, M, minT, histogram
h = setup_subplots(3,2,'vertical');
isub = 1;

hca = h(isub); isub = isub + 1;
nData = numel(out.M_bicmin);
histogram(hca,out.M_bicmin,0.5:1:8)
hold(hca,'on')
histogram(hca,out.M_bicchord,0.5:1:8)
hold(hca,'off')
hca.XLabel.String = 'Number of components';
hca.YLabel.String = 'Counts';
hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:200:5000;
%eT = logspace(1.5,4,10);
nData = numel(out.M_bicmin);
histogram(hca,out.minT_bicmin,eT,'Normalization','count')
hold(hca,'on')
histogram(hca,out.minT_bicchord,eT,'Normalization','count')
hold(hca,'off')
%hca.XScale = 'log';
hca.XLabel.String = 'min(T) (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.98 0.98])


hca = h(isub); isub = isub + 1;
eT = linspace(0,100,10000);
eT = 0:200:5000;
%eT = logspace(1.5,4,10);
nData = numel(out.M_bicmin);
histogram(hca,out.minT_bicmin,eT,'Normalization','cdf')
hold(hca,'on')
histogram(hca,out.minT_bicchord,eT,'Normalization','cdf')
hold(hca,'off')
%hca.XScale = 'log';
hca.XLabel.String = 'min(T) (eV)';
hca.YLabel.String = 'CDF';
hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eW = 0:0.05:1;
nData = numel(out.M_bicmin);
boxchart(hca,out.M_bicmin-0.1,out.minT_bicmin)
hold(hca,'on')
boxchart(hca,out.M_bicchord+0.1,out.minT_bicchord)
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
hca.YScale = 'log'; hca.YLim(1) = 10;
hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';

if 1
hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eW = 0:0.05:1;
nData = numel(out.M_bicmin);
boxchart(hca,out.M_bicmin-0.1,out.nMinT_bicmin)
hold(hca,'on')
boxchart(hca,out.M_bicchord+0.1,out.nMinT_bicchord)
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'Weight @min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
%hca.YScale = 'log'; hca.YLim(1) = 0.01;
hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';
end

if 1 %  (n,T)@min(T), scatter
  hca = h(isub); isub = isub + 1;
  eT = linspace(0,100,4000);
  eW = 0:0.05:1;
  nData = numel(out.M_bicmin);  
  scatter(hca,out.nMinT_bicmin,out.minT_bicmin,'.')
  hold(hca,'on')
  scatter(hca,out.nMinT_bicchord,out.minT_bicchord,'.')
  hold(hca,'off')
  hca.XLabel.String = 'Weight @min(T)';
  hca.YLabel.String = 'min(T) (eV)';
  hca.XScale = 'log'; hca.XLim(1) = 0.001;
  hca.YScale = 'log'; hca.YLim(1) = 10;
  hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
  irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])
  hca.Box = 'on';
end

if 0 %  (n,T)@min(T), contour
  hca = h(isub); isub = isub + 1;
  eT = linspace(0,1e3,20);
  eW = 0.00:0.05:1;
  eT = logspace(1,4,20);
  eW = logspace(-3,0,21);
  
  nData = numel(out.M_bicmin);
  [N edges mid loc] = histcn([out.nMinT_bicmin,out.minT_bicmin], eW, eT);
  contour(hca,mid{1},mid{2},N',17,'b','LineWidth',1)
  hold(hca,'on')
  [N edges mid loc] = histcn([out.nMinT_bicchord,out.minT_bicchord], eW, eT);
  contour(hca,mid{1},mid{2},N',17,'r','LineWidth',1)
  %scatter(hca,out.nMinT_bicchord,out.minT_bicchord,'.')
  hold(hca,'off')
  hca.XLabel.String = 'Weight @min(T)';
  hca.YLabel.String = 'min(T) (eV)';
  hca.XScale = 'log'; %hca.XLim(1) = 0.001;
  hca.YScale = 'log'; %hca.YLim(1) = 10;
  hca.XLim = eW([1 end]);
  hca.YLim = eT([1 end]);
  hca.Title.String = "nEvents = " + out.nEvents + ", nVDFs = " + out.nVDFs;
  irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])
  hca.Box = 'on';
end