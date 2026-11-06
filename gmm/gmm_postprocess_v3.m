% Save path
localuser = 'cecilianorgren';
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v1/'];
listFiles = dir([pathSave '*iFF*']);
nFiles = numel(listFiles);


tic
%stats = gmmCrossEventStatistics(listFiles(1:100),'stats',{'bicmin','bicchord','l1'},'thrL1',[0.05:0.05:0.2]);
%stats = gmmCrossEventStatistics(listFiles(1:2),'stats',{'bicmin','bicchord','l1','cold'},'thrL1',[0.05:0.05:0.2],'coldArgs',{'maxColdEV',500});
stats = gmmCrossEventStatistics(listFiles(1:20),'stats',{'bicmin','bicchord','l1'},'thrL1',[0.05:0.05:0.2],'coldArgs',{'maxColdEV',200});
toc
% perEventMeans = accumarray(out.bicmin.pre.event,out.bicmin.pre.minT,[],@mean)

%% Figure, M, minT, histogram, comparing pre post and all
h = setup_subplots(3,2,'vertical');
isub = 1;


method = 'l1';
iL1 = 1;

perEvent.pre.nVDFs = accumarray(stats.(method).pre.event,1,[],@sum);
perEvent.post.nVDFs = accumarray(stats.(method).post.event,1,[],@sum);
perEvent.all.nVDFs = accumarray(stats.(method).all.event,1,[],@sum);


nDataPre = numel(stats.(method).pre.M);
nDataPost = numel(stats.(method).post.M);
nDataAll = numel(stats.(method).all.M);

hca = h(isub); isub = isub + 1;
N = struct('pre',[],'post',[],'all',[]);
eT = 0.5:1:8;
cT = eT(1:end-1)+0.5*diff(eT);
histNorm = 'pdf';
N.pre = histcounts(stats.(method).pre.M(:,iL1),eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.M(:,iL1),eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.M(:,iL1),eT,'Normalization',histNorm);
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
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.98])


hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:100:3000;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
histNorm = 'count';
N.pre = histcounts(stats.(method).pre.minT(:,iL1),eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.minT(:,iL1),eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.minT(:,iL1),eT,'Normalization',histNorm);
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

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:100:3000;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
N.pre = histcounts(stats.(method).pre.minT(:,iL1),eT,'Normalization','cdf');
N.post = histcounts(stats.(method).post.minT(:,iL1),eT,'Normalization','cdf');
N.all = histcounts(stats.(method).all.minT(:,iL1),eT,'Normalization','cdf');
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
hca.YLim = [0 1];

hca = h(isub); isub = isub + 1; % delete(hca)
eT = linspace(0,100,4000);
eT = 0:100:2000;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
histNorm = 'count';
N.pre = histcounts(perEvent.pre.nVDFs,eT,'Normalization',histNorm);
N.post = histcounts(perEvent.post.nVDFs,eT,'Normalization',histNorm);
N.all = histcounts(perEvent.all.nVDFs,eT,'Normalization',histNorm);
steps(hca,eT,N.pre,'LineWidth',1)
hold(hca,'on')
steps(hca,eT,N.post,'LineWidth',1)
steps(hca,eT,N.all,'LineWidth',1)
hold(hca,'off')
%plot(hca,cT,N.pre,cT,N.post)
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XMinorGrid = 'on';
hca.XLabel.String = 'nVDFs/event';
hca.YLabel.String = histNorm;
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.98 0.98])



hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:0.05:1;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
histNorm = 'count';
N.pre = histcounts(stats.(method).pre.wMinT(:,iL1),eT,'Normalization',histNorm);
N.post = histcounts(stats.(method).post.wMinT(:,iL1),eT,'Normalization',histNorm);
N.all = histcounts(stats.(method).all.wMinT(:,iL1),eT,'Normalization',histNorm);
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
hca.YLabel.String = histNorm;
irf_legend(hca,{"pre-DF (N=" + stats.nVDFs.pre + ")","post-DF (N=" + stats.nVDFs.post + ")","all (N=" + stats.nVDFs.all + ")"}',[0.9 0.98])

hca = h(isub); isub = isub + 1;
eT = linspace(0,100,4000);
eT = 0:0.05:1;
cT = eT(1:end-1)+0.5*diff(eT);
N = struct('pre',[],'post',[],'all',[]);
N.pre = histcounts(stats.(method).pre.wMinT(:,iL1),eT,'Normalization','cdf');
N.post = histcounts(stats.(method).post.wMinT(:,iL1),eT,'Normalization','cdf');
N.all = histcounts(stats.(method).all.wMinT(:,iL1),eT,'Normalization','cdf');
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
hca.YLim = [0 1];

if strcmp(method,'l1')
  h(1).Title.String = {"nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.all,method + " < " + stats.thrL1(iL1)};
else
  h(1).Title.String = {"nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.all,method};
end

%% Figure, comparing bicmin, bicchord, L1.  --  pre post and all separately
h = setup_subplots(2,3,'vertical');
isub = 1;


method = 'l1';
iL1 = 1;

perEvent.pre.M = accumarray(stats.(method).pre.event,1,[],@sum);
perEvent.post.M = accumarray(stats.(method).post.event,1,[],@sum);
perEvent.all.M = accumarray(stats.(method).all.event,1,[],@sum);


nDataPre = numel(stats.(method).pre.M);
nDataPost = numel(stats.(method).post.M);
nDataAll = numel(stats.(method).all.M);

marker = 'none';

for where = ["pre","post","all"]  
  hca = h(isub); isub = isub + 1;
  N = struct('bicmin',[],'bicchord',[],'l1',[]);
  eT = 0.5:1:8;
  cT = eT(1:end-1)+0.5*diff(eT);
  histNorm = 'pdf';
  N.bicmin.(where) = histcounts(stats.('bicmin').(where).M,eT,'Normalization',histNorm);
  N.bicchord.(where) = histcounts(stats.('bicchord').(where).M,eT,'Normalization',histNorm);
  for iL1 = 1:numel(stats.thrL1)
    N.l1.(where)(:,iL1) = histcounts(stats.('l1').(where).M(:,iL1),eT,'Normalization',histNorm);
  end
  if 0 % plot as steps
    steps(hca,eT,N.bicmin.(where),'LineWidth',1)
    hold(hca,'on')
    steps(hca,eT,N.bicchord.(where),'LineWidth',1)
    for iL1 = 1:numel(stats.thrL1)
      steps(hca,eT,N.l1.(where)(:,iL1),'LineWidth',1)
    end
  else % plot line  
    plot(hca,cT,N.bicmin.(where),'LineWidth',1,'Marker',marker)
    hold(hca,'on')
    plot(hca,cT,N.bicchord.(where),'LineWidth',1,'Marker',marker)
    for iL1 = 1:numel(stats.thrL1)
      plot(hca,cT,N.l1.(where)(:,iL1),'LineWidth',1,'Marker',marker)
    end
  end
  hold(hca,'off')
  %plot(hca,0.5:1:8,N.pre,1:7,N.post)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'M';
  hca.XLabel.String = 'Number of components';
  hca.YLabel.String = histNorm;
  hca.XTick = 0:20;
  hca.XLim = eT([1 end]);
  irf_legend(hca,cat(2,"bicmin","bicchord","L1<" + stats.thrL1)',[1.02 0.98])
  hca.Title.String = {"nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(where) + ", where = " + where};
  
  
  hca = h(isub); isub = isub + 1;
  N = struct('bicmin',[],'bicchord',[],'l1',[]);
  eT = 0:100:3000;
  cT = eT(1:end-1)+0.5*diff(eT);
  histNorm = 'cdf';
  N.bicmin.(where) = histcounts(stats.('bicmin').(where).minT,eT,'Normalization',histNorm);
  N.bicchord.(where) = histcounts(stats.('bicchord').(where).minT,eT,'Normalization',histNorm);
  for iL1 = 1:numel(stats.thrL1)
    N.l1.(where)(:,iL1) = histcounts(stats.('l1').(where).minT(:,iL1),eT,'Normalization',histNorm);
  end
  steps(hca,eT,N.bicmin.(where),'LineWidth',1)
  hold(hca,'on')
  steps(hca,eT,N.bicchord.(where),'LineWidth',1)
  for iL1 = 1:numel(stats.thrL1)
    steps(hca,eT,N.l1.(where)(:,iL1),'LineWidth',1)
  end
  hold(hca,'off')
  %plot(hca,0.5:1:8,N.pre,1:7,N.post)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.XLabel.String = 'T_{min} (eV)';
  hca.XLabel.String = 'Temperature of coldest component (eV)';
  hca.YLabel.String = histNorm;
  %hca.XTick = 0:20;
  hca.XTick = 0:400:eT(end);
  hca.XLim = eT([1 end]);
  irf_legend(hca,cat(2,"bicmin","bicchord","L1<" + stats.thrL1)',[1.02 0.98])
  hca.XTickLabelRotation = 0;
  %hca.XLim = [0 2000];
  hca.YTick = 0:0.2:1;
  hca.YLim = [0 1];
  %hca.XMinorGrid = "on";
  % labels = hca.XTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.XTickLabels = labels;
  % labels = hca.YTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.YTickLabels = labels;
end
compact_panels(h,0.1,0.08,1); drawnow; compact_panels(h,0.1,0.08,1)

%% Figure, nColdComp  --  pre post and all separately
h = setup_subplots(2,3,'vertical');
isub = 1;


method = 'l1';
iL1 = 1;

perEvent.pre.M = accumarray(stats.(method).pre.event,1,[],@sum);
perEvent.post.M = accumarray(stats.(method).post.event,1,[],@sum);
perEvent.all.M = accumarray(stats.(method).all.event,1,[],@sum);


nDataPre = numel(stats.(method).pre.M);
nDataPost = numel(stats.(method).post.M);
nDataAll = numel(stats.(method).all.M);

marker = 'none';

for where = ["pre","post","all"]  
  hca = h(isub); isub = isub + 1;
  N = struct('bicmin',[],'bicchord',[],'l1',[]);
  eT = 0.5:1:8;
  cT = eT(1:end-1)+0.5*diff(eT);
  histNorm = 'pdf';
  N.bicmin.(where) = histcounts(stats.('bicmin').(where).M,eT,'Normalization',histNorm);
  N.bicchord.(where) = histcounts(stats.('bicchord').(where).M,eT,'Normalization',histNorm);
  for iL1 = 1:numel(stats.thrL1)
    N.l1.(where)(:,iL1) = histcounts(stats.('l1').(where).M(:,iL1),eT,'Normalization',histNorm);
  end
  if 0 % plot as steps
    steps(hca,eT,N.bicmin.(where),'LineWidth',1)
    hold(hca,'on')
    steps(hca,eT,N.bicchord.(where),'LineWidth',1)
    for iL1 = 1:numel(stats.thrL1)
      steps(hca,eT,N.l1.(where)(:,iL1),'LineWidth',1)
    end
  else % plot line  
    plot(hca,cT,N.bicmin.(where),'LineWidth',1,'Marker',marker)
    hold(hca,'on')
    plot(hca,cT,N.bicchord.(where),'LineWidth',1,'Marker',marker)
    for iL1 = 1:numel(stats.thrL1)
      plot(hca,cT,N.l1.(where)(:,iL1),'LineWidth',1,'Marker',marker)
    end
  end
  hold(hca,'off')
  %plot(hca,0.5:1:8,N.pre,1:7,N.post)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'M';
  hca.XLabel.String = 'Number of components';
  hca.YLabel.String = histNorm;
  hca.XTick = 0:20;
  hca.XLim = eT([1 end]);
  irf_legend(hca,cat(2,"bicmin","bicchord","L1<" + stats.thrL1)',[1.02 0.98])
  hca.Title.String = {"nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(where) + ", where = " + where};
  
  
  hca = h(isub); isub = isub + 1;
  N = struct('bicmin',[],'bicchord',[],'l1',[]);
  eT = -0.5:1:3;
  cT = eT(1:end-1)+0.5*diff(eT);
  histNorm = 'pdf';
  N.bicmin.(where) = histcounts(stats.('bicmin').(where).nCold,eT,'Normalization',histNorm);
  N.bicchord.(where) = histcounts(stats.('bicchord').(where).nCold,eT,'Normalization',histNorm);
  for iL1 = 1:numel(stats.thrL1)
    N.l1.(where)(:,iL1) = histcounts(stats.('l1').(where).nCold(:,iL1),eT,'Normalization',histNorm);
  end
  steps(hca,eT,N.bicmin.(where),'LineWidth',1)
  hold(hca,'on')
  steps(hca,eT,N.bicchord.(where),'LineWidth',1)
  for iL1 = 1:numel(stats.thrL1)
    steps(hca,eT,N.l1.(where)(:,iL1),'LineWidth',1)
  end
  hold(hca,'off')
  %plot(hca,0.5:1:8,N.pre,1:7,N.post)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.XLabel.String = 'T_{min} (eV)';
  hca.XLabel.String = 'nColdComp';
  hca.YLabel.String = histNorm;
  %hca.XTick = 0:20;
  %hca.XTick = 0:400:eT(end);
  hca.XLim = eT([1 end]);
  irf_legend(hca,cat(2,"bicmin","bicchord","L1<" + stats.thrL1)',[1.02 0.98])
  hca.XTickLabelRotation = 0;
  %hca.XLim = [0 2000];
  hca.YTick = 0:0.2:1;
  hca.YLim = [0 1];
  %hca.XMinorGrid = "on";
  % labels = hca.XTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.XTickLabels = labels;
  % labels = hca.YTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.YTickLabels = labels;
end
compact_panels(h,0.1,0.08,1); drawnow; compact_panels(h,0.1,0.08,1)

%% Figure, comparing bicmin, bicchord, L1.  --  pre post and all separately

 %sdf
method = 'l1';
iL1 = 1;
clear perEvent;
perEvent.pre.M = accumarray(stats.(method).pre.event,1,[],@sum);
perEvent.post.M = accumarray(stats.(method).post.event,1,[],@sum);
perEvent.all.M = accumarray(stats.(method).all.event,1,[],@sum);
%prctile25 = @(x)prctile(x,25);
for stat = ["bicmin","bicchord","l1"]
  for where = ["pre","post","all"]
    for quantity = ["minT"]
      for fun = {@min,@median,@mean,@max,@prctile10,@prctile25}
        switch stat
          case 'l1'
            for iL1 = 1:numel(stats.thrL1)
              arr(:,iL1) = accumarray(stats.(stat).(where).event,stats.(stat).(where).(quantity)(:,iL1),[],fun{:});
            end
            perEvent.(stat).(where).(char(fun{:})).(quantity) = arr;
          otherwise
            perEvent.(stat).(where).(char(fun{:})).(quantity) = accumarray(stats.(stat).(where).event,stats.(stat).(where).(quantity),[],fun{:});
        end
      end
    end
  end
end


nDataPre = numel(stats.(method).pre.M);
nDataPost = numel(stats.(method).post.M);
nDataAll = numel(stats.(method).all.M);


h = setup_subplots(3,4,'vertical');
isub = 1;

marker = 'none';

fun = "median";
for fun = ["median","mean","prctile10","prctile25"]  
for where = ["pre","post","all"]  
  
  hca = h(isub); isub = isub + 1;
  N = struct('bicmin',[],'bicchord',[],'l1',[]);
  eT = 0:100:3000;
  cT = eT(1:end-1)+0.5*diff(eT);
  histNorm = 'cdf';
  N.bicmin.(where) = histcounts(perEvent.('bicmin').(where).(fun).minT,eT,'Normalization',histNorm);
  N.bicchord.(where) = histcounts(perEvent.('bicchord').(where).(fun).minT,eT,'Normalization',histNorm);
  for iL1 = 1:numel(stats.thrL1)
    N.l1.(where)(:,iL1) = histcounts(perEvent.('l1').(where).(fun).minT(:,iL1),eT,'Normalization',histNorm);
  end
  steps(hca,eT,N.bicmin.(where),'LineWidth',1)
  hold(hca,'on')
  steps(hca,eT,N.bicchord.(where),'LineWidth',1)
  for iL1 = 1:numel(stats.thrL1)
    steps(hca,eT,N.l1.(where)(:,iL1),'LineWidth',1)
  end
  hold(hca,'off')
  %plot(hca,0.5:1:8,N.pre,1:7,N.post)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Title.String = "where = " + where;
  %hca.XLabel.String = 'T_{min} (eV)';
  %hca.XLabel.String = sprintf('%s temperature of coldest component per event (eV)',char(fun));
  hca.XLabel.String = sprintf('Per-event %s T_{coldest} (eV)',char(fun));
  hca.YLabel.String = histNorm;
  %hca.XTick = 0:20;
  hca.XTick = 0:400:eT(end);
  hca.XLim = eT([1 end]);
  irf_legend(hca,cat(2,"bicmin","bicchord","L1<" + stats.thrL1)',[1.02 0.98])
  hca.XTickLabelRotation = 0;
  %hca.XLim = [0 2000];
  hca.YTick = 0:0.2:1;
  hca.YLim = [0 1];
  %hca.XMinorGrid = "on";
  % labels = hca.XTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.XTickLabels = labels;
  % labels = hca.YTickLabels;
  % emptylbs = cellfun(@(x) '',labels,'UniformOutput',false);
  % labels(2:2:end) = emptylbs(2:2:end);
  % hca.YTickLabels = labels;
end
end
%compact_panels(h,0.1,0.08,1); drawnow; compact_panels(h,0.1,0.08,1)















%% Figure, M, minT, histogram, comparing bicmin and bicchord
h = setup_subplots(3,2,'vertical');
isub = 1;

region = 'post';

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
irf_legend(hca,{'BIC: min','BIC: chord'}',[0.02 0.98])


hca = h(isub); isub = isub + 1;
%eT = linspace(0,100,4000);
%eW = 0:0.05:1;
%nData = numel(stats.M_bicmin);
boxchart(hca,stats.('bicmin').(region).M(:,iL1)-0.15,stats.('bicmin').(region).minT(:,iL1))
hold(hca,'on')
boxchart(hca,stats.('bicchord').(region).M(:,iL1)+0.15,stats.('bicchord').(region).minT(:,iL1))
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
hca.YScale = 'log'; hca.YLim(1) = 10;
%hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';

if 1
hca = h(isub); isub = isub + 1;
%eT = linspace(0,100,4000);
%eW = 0:0.05:1;
%nData = numel(stats.M_bicmin);
boxchart(hca,stats.('bicmin').(region).M(:,iL1)-0.15,stats.('bicmin').(region).wMinT(:,iL1))
hold(hca,'on')
boxchart(hca,stats.('bicchord').(region).M(:,iL1)+0.15,stats.('bicchord').(region).wMinT(:,iL1))
hold(hca,'off')
hca.XLabel.String = 'M';
hca.YLabel.String = 'Weight @min(T) (eV)';
%hca.XScale = 'log'; hca.XLim(1) = 0.001;
%hca.YScale = 'log'; hca.YLim(1) = 0.01;
%hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
irf_legend(hca,{'BIC: min','BIC: chord'}',[1.02 0.98])
hca.Box = 'on';
end

if 1 %  (n,T)@min(T), scatter
  hca = h(isub); isub = isub + 1;
  %eT = linspace(0,100,4000);
  %eW = 0:0.05:1;
  %nData = numel(stats.M_bicmin);  
  scatter(hca,stats.('bicmin').(region).wMinT(:,iL1),stats.('bicmin').(region).minT(:,iL1),'.')
  hold(hca,'on')
  scatter(hca,stats.('bicchord').(region).wMinT(:,iL1),stats.('bicchord').(region).minT(:,iL1),'.')
  hold(hca,'off')
  hca.XLabel.String = 'Weight @min(T)';
  hca.YLabel.String = 'min(T) (eV)';
  hca.XScale = 'log'; hca.XLim(1) = 0.001;
  hca.YScale = 'log'; hca.YLim(1) = 10;
  %hca.Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs;
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


h(1).Title.String = "nEvents = " + stats.nEvents + ", nVDFs = " + stats.nVDFs.(region) + ", region = " + region;