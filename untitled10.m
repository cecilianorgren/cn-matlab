% plot omni data
%% Multiple events
tint_events = {irf.tint('2015-09-14T16:06:51.000/2015-09-14T16:07:33.000'),...
               irf.tint('2015-09-23T08:15:03.000/2015-09-23T08:15:36.000'),...
               irf.tint('2015-09-23T10:56:28.000/2015-09-23T10:56:58.000'),...
               irf.tint('2015-09-25T09:24:22.793/2015-09-25T09:24:52.977'),...
               irf.tint('2015-09-25T09:57:43.000/2015-09-25T09:58:28.000'),...
               irf.tint('2015-10-03T10:45:51.000/2015-10-03T10:46:34.000'),...
               irf.tint('2015-10-03T13:27:12.396/2015-10-03T13:27:42.044'),...
               irf.tint('2015-10-06T17:28:07.567/2015-10-06T17:28:17.015'),...
               irf.tint('2015-10-08T07:41:32.693/2015-10-08T07:41:59.549'),...
               irf.tint('2015-10-11T11:05:13.000/2015-10-11T11:05:57.000'),...
               irf.tint('2015-10-20T06:16:03.000/2015-10-20T06:17:26.000'),...
               irf.tint('2015-10-22T13:40:19.000/2015-10-22T13:40:48.000'),...
               irf.tint('2015-10-31T05:45:52.000/2015-10-31T05:46:27.000'),...
               irf.tint('2015-10-31T05:45:52.000/2015-10-31T05:46:27.000'),...
               irf.tint('2015-11-05T14:07:07.131/2015-11-05T14:07:44.639'),...
               irf.tint('2015-11-05T14:36:39.263/2015-11-05T14:36:44.695'),...
               irf.tint('2015-11-06T06:57:42.000/2015-11-06T06:58:26.000'),...
               irf.tint('2015-11-06T13:23:48.000/2015-11-06T13:24:29.000'),...
               irf.tint('2015-11-06T13:26:17.915/2015-11-06T13:26:32.647'),...
               irf.tint('2015-11-08T14:02:51.000/2015-11-08T14:03:23.000'),...
               irf.tint('2015-11-09T10:06:54.162/2015-11-09T10:07:02.160')};
nevents = numel(tint_events);

%h = irf_plot(2);
all_data = [];
for ievent = 1:nevents
  tint_event = tint_events{ievent};
  tint_omni = tint_event(1)+[-20*60 0];
  data = irf_get_data(tint_omni,'b,bx,bygsm,bzgsm,v,ts','omni_min');
  all_data = cat(3,all_data,data);
end

%%
By = squeeze(all_data(:,[4],:));
Bz = squeeze(all_data(:,[5],:));
v =  squeeze(all_data(:,[6],:));
ts =  squeeze(all_data(:,[7],:));
clock_angle = atand(Bz./By);

linestyles = {'-','--',':','.-'};

h = setup_subplots(4,1);
isub = 1;
hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
hca.LineStyleOrder = linestyles;
hold(hca,'on')
plot(hca,[-19:1:0],By)
hold(hca,'off')
hca.XLabel.String = 't-t_{start}';
hca.YLabel.String = 'B_y (nT)';
%irf_legend(hca,{'|B|','B_x','B_y','B_z'},[0.98 0.98])

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
hca.LineStyleOrder = linestyles;
hold(hca,'on')
plot(hca,[-19:1:0],Bz)
hold(hca,'off')
%plot(hca,squeeze(all_data(:,[4],:)))
hca.XLabel.String = 't-t_{start}';
hca.YLabel.String = 'B_z (nT)';

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
hca.LineStyleOrder = linestyles;
hold(hca,'on')
plot(hca,[-19:1:0],clock_angle)
hold(hca,'off')
%plot(hca,squeeze(all_data(:,[4],:)))
hca.XLabel.String = 't-t_{start}';
hca.YLabel.String = 'clock angle (deg)';

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
hca.LineStyleOrder = linestyles;
hold(hca,'on')
plot(hca,[-19:1:0],v)
hold(hca,'off')
%plot(hca,squeeze(all_data(:,[4],:)))
hca.XLabel.String = 't-t_{start}';
hca.YLabel.String = 'v (km/s)';

if 0 % ts
  hca = h(isub); isub = isub + 1;
  hca.ColorOrder = pic_colors('matlab');
  hca.LineStyleOrder = linestyles;
  hold(hca,'on')
  plot(hca,[-19:1:0],ts)
  hold(hca,'off')
  %plot(hca,squeeze(all_data(:,[4],:)))
  hca.XLabel.String = 't-t_{start}';
  hca.YLabel.String = 'time since/to bowshock (...)';
end

compact_panels(h,0.01)

for ip = 1:nrows*ncols
  h(ip).LineWidth = 1;
  %h(ip).LineStyleOrder = {'-','--','.-',':'};
end
%%
h = irf_plot(1);

hca = irf_panel('B');
irf_plot(hca,[all_data(:,1,1) squeeze(all_data(:,[4],:))])
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'|B|','B_x','B_y','B_z'},[0.98 0.98])

% hca = irf_panel('v');
% irf_plot(hca,data(:,[1 8]))
% hca.YLabel.String = 'v (km/s)';

%% Sigmle event
tint_event = irf.tint('2015-09-14T16:06:51.000/2015-09-14T16:07:33.000');
tint_omni = tint_event(1)+[-20*60 0];
data = irf_get_data(tint_omni,'b,bx,bygsm,bzgsm,ts','omni_min');

h = irf_plot(2);

hca = irf_panel('B');
irf_plot(hca,data(:,1:5))
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'|B|','B_x','B_y','B_z'},[0.98 0.98])

hca = irf_panel('v');
irf_plot(hca,data(:,[1 8]))
hca.YLabel.String = 'v (km/s)';
