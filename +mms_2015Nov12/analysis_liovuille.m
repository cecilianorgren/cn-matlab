% analysis_liovuille

% particles starting on left and exiting on right
leftstart = find(allz0<0); % index
leftstop = find(allzend<0);
rightstart = find(allz0>0);
rightstop = find(allzend>0);
centerpass = find(allzmin<10e3); 

leftleft = intersect(leftstart,leftstop);
leftright = intersect(leftstart,rightstop);
rightleft = intersect(rightstart,leftstop);
rightright = intersect(rightstart,rightstop);

leftcenterleft = intersect(leftleft,centerpass);
leftcenterright = intersect(leftright,centerpass);
rightcenterleft = intersect(rightleft,centerpass);
rightcenterright = intersect(rightright,centerpass);


leftleft_NE = histc(eENERGY(leftleft),edgesE);
leftright_NE = histc(eENERGY(leftright),edgesE);
rightleft_NE = histc(eENERGY(rightleft),edgesE);
rightright_NE = histc(eENERGY(rightright),edgesE);

leftcenterleft_NE = histc(eENERGY(leftcenterleft),edgesE);
leftcenterright_NE = histc(eENERGY(leftcenterright),edgesE);
rightcenterleft_NE = histc(eENERGY(rightcenterleft),edgesE);
rightcenterright_NE = histc(eENERGY(rightcenterright),edgesE);

limit = 2;
edgesDiffE = [-limit:0.2:limit];
diffE = (allEend-allE0)./allE0;
leftleft_NEdiff = histc(diffE(leftleft),edgesDiffE);
leftright_NEdiff = histc(diffE(leftright),edgesDiffE);
rightleft_NEdiff = histc(diffE(rightleft),edgesDiffE);
rightright_NEdiff = histc(diffE(rightright),edgesDiffE);

f_edges = logspace(-30,-25,20);
left_f0 = histc(f0(leftstart),f_edges);
right_f0 = histc(f0(rightstart),f_edges);

z_min_edges = (-30:10:30)*1e3;
left_zmin = histc(allzmin(leftstart),z_min_edges);
right_zmin = histc(allzmin(rightstart),z_min_edges);


figure(61)
nrows = 4;
ncols = 3;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;
 
if 1 % start and stop side
  hca = h(isub); isub = isub + 1;
  %hb = bar(hca,log10(edgesE),[leftleft_NE leftright_NE rightleft_NE rightright_NE],'stacked'); colormap('cool')
  hb = bar(hca,log10(edgesE),[leftcenterleft_NE leftcenterright_NE rightcenterleft_NE rightcenterright_NE],'stacked'); colormap(mms_colors('1234'))
  legend(hca,'leftleft','leftright','rightleft','rightright','location','eastoutside')
end
if 1 % start and stop energy
  hca = h(isub); isub = isub + 1;
  colorE = zeros(numel(allEend),3);
  colors = mms_colors('1234');
  colorE(leftleft,:) = repmat(colors(1,:),numel(leftleft),1);
  colorE(leftright,:) = repmat(colors(2,:),numel(leftright),1);
  colorE(rightleft,:) = repmat(colors(3,:),numel(rightleft),1);
  colorE(rightright,:) = repmat(colors(4,:),numel(rightright),1);
  hscat = scatter(hca,allE0,allEend,allE0*0+20,colorE);
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.XLabel.String = 'E_{start}';
  hca.YLabel.String = 'E_{end}';
  axis(hca,'square')
  hca.XLim = [1e0 1e4];
  hca.YLim = [1e0 1e4];
  hold(hca,'on')
  plot(hca,hca.XLim,hca.YLim,'k')
  text(hca,1e1,1e3,'E increase')
  text(hca,1e3,1e1,'E decrease')
  hold(hca,'off')
end
if 1 % E0 vs allTend
  hca = h(isub); isub = isub + 1;
  colorE = zeros(numel(allEend),3);
  colors = mms_colors('1234');
  colorE(leftleft,:) = repmat(colors(1,:),numel(leftleft),1);
  colorE(leftright,:) = repmat(colors(2,:),numel(leftright),1);
  colorE(rightleft,:) = repmat(colors(3,:),numel(rightleft),1);
  colorE(rightright,:) = repmat(colors(4,:),numel(rightright),1);
  hscat = scatter(hca,allTend,allE0,diffE*0+20,colorE);  
  hca.XLabel.String = 'T_{end}';
  hca.YLabel.String = 'E_{start}';
  axis(hca,'square')
  hca.XLim = [0 1];
  %hca.YLim = [-2 3];
end
if 1 % diff E vs allTend
  hca = h(isub); isub = isub + 1;
  colorE = zeros(numel(allEend),3);
  colors = mms_colors('1234');
  colorE(leftleft,:) = repmat(colors(1,:),numel(leftleft),1);
  colorE(leftright,:) = repmat(colors(2,:),numel(leftright),1);
  colorE(rightleft,:) = repmat(colors(3,:),numel(rightleft),1);
  colorE(rightright,:) = repmat(colors(4,:),numel(rightright),1);
  hscat = scatter(hca,allTend,diffE,diffE*0+20,colorE);  
  hca.XLabel.String = 'T_{end}';
  hca.YLabel.String = '(E_{end}-E_{start})/E_{start}';
  axis(hca,'square')
  hca.XLim = [0 1];
  %hca.YLim = [-2 3];
end
if 1 % diff E vs E0
  hca = h(isub); isub = isub + 1;
  colorE = zeros(numel(allEend),3);
  colors = mms_colors('1234');
  colorE(leftleft,:) = repmat(colors(1,:),numel(leftleft),1);
  colorE(leftright,:) = repmat(colors(2,:),numel(leftright),1);
  colorE(rightleft,:) = repmat(colors(3,:),numel(rightleft),1);
  colorE(rightright,:) = repmat(colors(4,:),numel(rightright),1);
  hscat = scatter(hca,allE0,diffE,diffE*0+20,colorE);  
  hca.XLabel.String = 'E_{0}';
  hca.YLabel.String = '(E_{end}-E_{start})/E_{start}';
  axis(hca,'square')
  %hca.XLim = [0 1];
  %hca.YLim = [-2 3];
end
if 1 % diff E
  hca = h(isub); isub = isub + 1;
  limit = 3;
  %hist(hca,diffE,edgesDiffE)
  hb = bar(hca,edgesDiffE,[leftleft_NEdiff leftright_NEdiff rightleft_NEdiff rightright_NEdiff],'stacked'); 
  hca.XLim = [-limit, limit];
  hca.XLabel.String = '(E_{end}-E_{start})/E_{start}';
  hca.YLabel.String = 'Counts';
  legend(hca,'leftleft','leftright','rightleft','rightright','location','eastoutside')
end
if 1 % integration time
  hca = h(isub); isub = isub + 1;
  hist(hca,allTend,[0:0.1:1.9])
  hca.XLabel.String = 'T_{int}';
  hca.YLabel.String = 'counts';
end
if 1 % f0
  hca = h(isub); isub = isub + 1;
  colors = mms_colors('matlab');
  hca.ColorOrder = colors;
  bar(hca,log10(f_edges),[left_f0 right_f0]); colormap(hca,colors(1:2,:));
  hca.XLim = log10(f_edges([1 end]));
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hca.ColorOrder = colors;
  irf_legend(hca,{'left','right'},[0.96 0.95])
  %legend(hca,'left','right','location','eastoutside');  
end
if 1 % zmin
  hca = h(isub); isub = isub + 1;
  colors = mms_colors('matlab');
  hca.ColorOrder = colors;
  bar(hca,z_min_edges,[left_zmin right_zmin]); colormap(hca,colors(1:2,:));
  %hca.XLim = log10(f_edges([1 end]));
  %hca.XScale = 'log';
  %hca.YScale = 'log';
  hca.ColorOrder = colors;
  irf_legend(hca,{'left','right'},[0.96 0.95])
  %legend(hca,'left','right','location','eastoutside');  
end

if 1 % pitchangle distributions of obs and map, left
  hca = h(isub); isub = isub + 1;
  
  colors = mms_colors('matlab'); colors = colors(1:3,:);
  time1 = irf_time('2015-11-12T07:19:20.900Z','utc>EpochTT');   
  time2 = t_left;
  tind1_obs  = find(abs(obsPDist.time-time1)==min(abs(obsPDist.time-time1)));
  tind1_map  = find(abs(tsFmap.time-time1)==min(abs(tsFmap.time-time1)));
  
  pitch1_map = tsFmap(tind1_map).pitchangles(gseB1,15);
  pitch1_obs = obsPDist(tind1_obs).pitchangles(gseB1,15);
  
  hca.ColorOrder = colors;
  plot(hca,pitch1_obs.depend{1}(1,:),squeeze(pitch1_obs.data(1,:,[1 ceil(numel(pitch1_obs.depend{2})/2) numel(pitch1_obs.depend{2})])));
  hold(hca,'on')
  hca.ColorOrder = colors;
  hca.LineStyleOrder = '--';
  plot(hca,pitch1_map.depend{1}(1,:),squeeze(pitch1_map.data(1,:,[1 ceil(numel(pitch1_map.depend{2})/2) numel(pitch1_map.depend{2})])));
  hold(hca,'off')
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' pitch1_obs.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 1000];
  legend(hca,{'0','90','180'})
  hca.YTick = 10.^[-3:5]*1e-30;
  hca.YLim = [1e-1 1e5]*1e-30;
  hca.Title.String = pitch1_obs.time.utc;
end
if 1 % pitchangle distributions of obs and map, center
  hca = h(isub); isub = isub + 1;
  
  colors = mms_colors('matlab'); colors = colors(1:3,:);
  time2 = irf_time('2015-11-12T07:19:21.500Z','utc>EpochTT'); 
  time2 = t_center;
  tind2_obs  = find(abs(obsPDist.time-time2)==min(abs(obsPDist.time-time2)));
  tind2_map  = find(abs(tsFmap.time-time2)==min(abs(tsFmap.time-time2)));
  
  pitch2_map = tsFmap(tind2_map).pitchangles(gseB1,15);
  pitch2_obs = obsPDist(tind2_obs).pitchangles(gseB1,15);
  
  hca.ColorOrder = colors;
  plot(hca,pitch2_obs.depend{1}(1,:),squeeze(pitch2_obs.data(1,:,[1 ceil(numel(pitch2_obs.depend{2})/2) numel(pitch2_obs.depend{2})])));
  hold(hca,'on')
  hca.ColorOrder = colors;
  hca.LineStyleOrder = '--';
  plot(hca,pitch2_map.depend{1}(1,:),squeeze(pitch2_map.data(1,:,[1 ceil(numel(pitch2_map.depend{2})/2) numel(pitch2_map.depend{2})])));
  hold(hca,'off')
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' pitch1_obs.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 1000];
  legend(hca,{'0','90','180'})
  hca.YTick = 10.^[-3:5]*1e-30;
  hca.YLim = [1e-1 1e5]*1e-30;
  hca.Title.String = pitch2_obs.time.utc;
end
if 1 % pitchangle distributions of obs and map, right
  hca = h(isub); isub = isub + 1;
  
  colors = mms_colors('matlab'); colors = colors(1:3,:);
  time2 = irf_time('2015-11-12T07:19:21.500Z','utc>EpochTT'); 
  time2 = t_right;
  tind2_obs  = find(abs(obsPDist.time-time2)==min(abs(obsPDist.time-time2)));
  tind2_map  = find(abs(tsFmap.time-time2)==min(abs(tsFmap.time-time2)));
  
  pitch2_map = tsFmap(tind2_map).pitchangles(gseB1,15);
  pitch2_obs = obsPDist(tind2_obs).pitchangles(gseB1,15);
  
  hca.ColorOrder = colors;
  plot(hca,pitch2_obs.depend{1}(1,:),squeeze(pitch2_obs.data(1,:,[1 ceil(numel(pitch2_obs.depend{2})/2) numel(pitch2_obs.depend{2})])));
  hold(hca,'on')
  hca.ColorOrder = colors;
  hca.LineStyleOrder = '--';
  plot(hca,pitch2_map.depend{1}(1,:),squeeze(pitch2_map.data(1,:,[1 ceil(numel(pitch2_map.depend{2})/2) numel(pitch2_map.depend{2})])));
  hold(hca,'off')
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' pitch1_obs.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 1000];
  legend(hca,{'0','90','180'})
  hca.YTick = 10.^[-3:5]*1e-30;
  hca.YLim = [1e-1 1e5]*1e-30;
  hca.Title.String = pitch2_obs.time.utc;
end