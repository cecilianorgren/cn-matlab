% Save path
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v1/'];
listFiles = dir([pathSave '*iFF*']);
nFiles = numel(listFiles);


%out = gmmCrossEventStatistics(listFiles(1:50),'version',{'minbick','chordbick'});
stats = gmmCrossEventStatistics(listFiles(1:50));
stats2 = gmmCrossEventStatistics(listFiles(1:10),'stats',{'minbick','chordbick','l1'},'thrL1',[0.1 0.2]);

% perEventMeans = accumarray(out.bicmin.pre.event,out.bicmin.pre.minT,[],@mean)
%% Figure, M, minT, histogram, comparing pre post and all
h = setup_subplots(3,2,'vertical');
isub = 1;

nDataPre = numel(stats.(method).pre.M);
nDataPost = numel(stats.(method).post.M);
nDataAll = numel(stats.(method).all.M);


method = 'bicchord';
%method = 'bicmin';
hca = h(isub); isub = isub + 1;
N = struct('pre',[],'post',[],'all',[]);
eT = 0.5:1:8;
cT = eT(1:end-1)+0.5*diff(eT);
histNorm = 'count';
N.pre = histcounts(stats.(method).pre.M,eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.M,eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.M,eT,'Normalization',histNorm);
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,0.5:1:8,N.pre,1:7,N.post)
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'M';
hca.YLabel.String = histNorm;
hca.XTick = 0:20;
hca.XLim = eT([1 end]);
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.4 0.98])


%method = 'bicmin';
hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:100:3000;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
histNorm = 'count';
N.pre = histcounts(stats.(method).pre.minT,eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.minT,eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.minT,eT,'Normalization',histNorm);
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,cT,N.pre,cT,N.post)
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XMinorGrid = 'on';
hca.XLabel.String = 'T_{min} (eV)';
hca.YLabel.String = histNorm;
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.98])

%method = 'bic';
hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:100:3000;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
N.pre = histcounts(stats.(method).pre.minT,eT,'Normalization','cdf');
N.post = histcounts(stats.(method).post.minT,eT,'Normalization','cdf');
N.all = histcounts(stats.(method).all.minT,eT,'Normalization','cdf');
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,cT,N.pre,cT,N.post)
hca.YLim(1) = 0;
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XMinorGrid = 'on';
%hca.XTick = 0:200:eT(end);
hca.XMinorTick = 'on';% 0:200:eT(end);
hca.XLabel.String = 'T_{min} (eV)';
hca.YLabel.String = 'CDF';
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.1])

hca = h(isub); isub = isub + 1;
delete(hca)

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:0.05:1;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
histNorm = 'count';
N.pre = histcounts(stats.(method).pre.wMinT,eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.wMinT,eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.wMinT,eT,'Normalization',histNorm);
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,cT,N.pre,cT,N.post)
hca.YLim(1) = 0;
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XMinorGrid = 'on';
%hca.XTick = 0:200:eT(end);
hca.XMinorTick = 'on';% 0:200:eT(end);
hca.XLabel.String = 'weight @ T_{min}';
hca.YLabel.String = 'CDF';
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.98])

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:0.05:1;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
N.pre = histcounts(stats.(method).pre.wMinT,eT,'Normalization','cdf');
N.post = histcounts(stats.(method).post.wMinT,eT,'Normalization','cdf');
N.all = histcounts(stats.(method).all.wMinT,eT,'Normalization','cdf');
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,cT,N.pre,cT,N.post)
hca.YLim(1) = 0;
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XMinorGrid = 'on';
%hca.XTick = 0:200:eT(end);
hca.XMinorTick = 'on';% 0:200:eT(end);
hca.XLabel.String = 'weight @ T_{min}';
hca.YLabel.String = 'CDF';
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.1])

h(1).Title.String = {"nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.all,method};

%%

%% Figure, M, minT, histogram
h = setup_subplots(3,2,'vertical');
isub = 1;

region = 'pre';
hca = h(isub); isub = isub + 1;
nData = numel(stats.bicmin.(region).M);
histogram(hca,stats.bicmin.(region).M,0.5:1:8,'Normalization','pdf')
hold(hca,'on')
histogram(hca,stats.bicchord.(region).M,0.5:1:8,'Normalization','pdf')
hold(hca,'off')
hca.XLabel.String = 'Number of components';
hca.YLabel.String = 'Counts';
hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(region);
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])


hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:200:5000;
%eT = logspace(1.5,4,10);
nData = numel(stats.bicmin.(region).minT);
histogram(hca,stats.bicmin.(region).minT,eT,'Normalization','pdf')
hold(hca,'on')
histogram(hca,stats.bicchord.(region).minT,eT,'Normalization','pdf')
hold(hca,'off')
%hca.XScale = 'log';
hca.XLabel.String = 'min(T) (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(region);
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.98 0.98])

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:200:5000;
%eT = logspace(1.5,4,10);
nData = numel(stats.bicmin.(region).minT);
histogram(hca,stats.bicmin.(region).minT,eT,'Normalization','cdf')
hold(hca,'on')
histogram(hca,stats.bicchord.(region).minT,eT,'Normalization','cdf')
hold(hca,'off')
%hca.XScale = 'log';
hca.XLabel.String = 'min(T) (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(region);
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])

%%
hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eW = 0:0.05:1;
nData = numel(stats.M_bicmin);
boxchart(hca,stats.M_bicmin-0.1,stats.minT_bicmin)
hold(hca,'on')
boxchart(hca,stats.M_bicchord+0.1,stats.minT_bicchord)
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
hca.YScale = 'log'; hca.YLim(1) = 10;
hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';

if 1
hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eW = 0:0.05:1;
nData = numel(stats.M_bicmin);
boxchart(hca,stats.M_bicmin-0.1,stats.nMinT_bicmin)
hold(hca,'on')
boxchart(hca,stats.M_bicchord+0.1,stats.nMinT_bicchord)
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'Weight @min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
%hca.YScale = 'log'; hca.YLim(1) = 0.01;
hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';
end

if 1 %  (n,T)@min(T), scatter
  hca = h(isub); isub = isub + 1;
  eT = linspace(0,100,4000);
  eW = 0:0.05:1;
  nData = numel(stats.M_bicmin);  
  scatter(hca,stats.nMinT_bicmin,stats.minT_bicmin,'.')
  hold(hca,'on')
  scatter(hca,stats.nMinT_bicchord,stats.minT_bicchord,'.')
  hold(hca,'off')
  hca.XLabel.String = 'Weight @min(T)';
  hca.YLabel.String = 'min(T) (eV)';
  hca.XScale = 'log'; hca.XLim(1) = 0.001;
  hca.YScale = 'log'; hca.YLim(1) = 10;
  hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
  irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])
  hca.Box = 'on';
end

if 0 %  (n,T)@min(T), contour
  hca = h(isub); isub = isub + 1;
  eT = linspace(0,1e3,20);
  eW = 0.00:0.05:1;
  eT = logspace(1,4,20);
  eW = logspace(-3,0,21);
  
  nData = numel(stats.M_bicmin);
  [N edges mid loc] = histcn([stats.nMinT_bicmin,stats.minT_bicmin], eW, eT);
  contour(hca,mid{1},mid{2},N',17,'b','LineWidth',1)
  hold(hca,'on')
  [N edges mid loc] = histcn([stats.nMinT_bicchord,stats.minT_bicchord], eW, eT);
  contour(hca,mid{1},mid{2},N',17,'r','LineWidth',1)
  %scatter(hca,out.nMinT_bicchord,out.minT_bicchord,'.')
  hold(hca,'off')
  hca.XLabel.String = 'Weight @min(T)';
  hca.YLabel.String = 'min(T) (eV)';
  hca.XScale = 'log'; %hca.XLim(1) = 0.001;
  hca.YScale = 'log'; %hca.YLim(1) = 10;
  hca.XLim = eW([1 end]);
  hca.YLim = eT([1 end]);
  hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
  irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])
  hca.Box = 'on';
end