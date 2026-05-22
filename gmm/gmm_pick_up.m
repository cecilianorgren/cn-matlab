

c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('tsVExB = gseVExB?.resample(tsV{1});',ic)

figure(77)
h = setup_subplots(3,1);

isub = 1;
hca = h(isub); isub = isub + 1;

iK = 1;
vx = tsV{iK}.x.data;
vy = tsV{iK}.y.data;
ss = tsT{iK}.trace.data/3*0+5;
scatter(hca,vx,vy,ss)
hold(hca,'on')
iK = 2;
vx = tsV{iK}.x.data;
vy = tsV{iK}.y.data;
ss = tsT{iK}.trace.data/3*0+5;
scatter(hca,vx,vy,ss)
hold(hca,'off')
axis(hca,'equal')
axis(hca,'square')


hca = h(isub); isub = isub + 1;

iK = 1;
vx = tsV{iK}.x.data - tsVExB.x.data;
vy = tsV{iK}.y.data - tsVExB.y.data;
ss = tsT{iK}.trace.data/3*0+5;
scatter(hca,vx,vy,ss)
hold(hca,'on')
iK = 2;
vx = tsV{iK}.x.data - tsVExB.x.data;
vy = tsV{iK}.y.data - tsVExB.y.data;
ss = tsT{iK}.trace.data/3*0+5;
scatter(hca,vx,vy,ss)
hold(hca,'off')
axis(hca,'equal')
axis(hca,'square')



hca = h(isub); isub = isub + 1;

iK = 1;
vx1 = tsV{iK}.x.data - tsVExB.x.data;
vy1 = tsV{iK}.y.data - tsVExB.y.data;

iK = 2;
vx2 = tsV{iK}.x.data - tsVExB.x.data;
vy2 = tsV{iK}.y.data - tsVExB.y.data;
histcn_plot([vx1,vy1;vx2,vy2],[-1000:100:1000],[-1000:100:1000])
axis(hca,'equal')
axis(hca,'square')


c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))