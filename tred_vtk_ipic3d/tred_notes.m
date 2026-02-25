%% tred_find_cores
core_satellite_position=tred_find_cores([5 9],[6 8],[6 9]);
limits=[min(core_satellite_position{3}(:,1)) max(core_satellite_position{3}(:,1));...
       min(core_satellite_position{3}(:,2)) max(core_satellite_position{3}(:,2));...
       min(core_satellite_position{3}(:,3)) max(core_satellite_position{3}(:,3))]

%% tred_plot_pos
tred_plot_pos(core_satellite_position{1},core_satellite_position{2},core_satellite_position{3})

%% tred_plot_cores (does same as tred_plot_pos)
core_satellite_position=tred_find_cores([5 9],[6 8],[6 9]);
figure;tred_plot_cores(core_satellite_position{1},core_satellite_position{2},core_satellite_position{3})

%% tred_load_n_plot_ov
%cores_to_plot=unique(core_satellite_position{1});
%ind=find(interp_sat{4}<0.5);
%cores_to_plot=interp_sat{1}(ind,:);
%sat_to_plot=interp_sat{2}(ind,:);

cores_to_plot=unique((core_satellite_position{1}(:)));
sat_to_plot=(core_satellite_position{2}(:));

for k=1:length(cores_to_plot)
c=cores_to_plot(k);s=14;%sat_to_plot(k);
figure(33);tred_load_n_plot_ov(c,s);
userdata=get(gcf,'userdata');
irf_zoom(userdata.subplot_handles,'x',[25 45]);
set(gcf,'PaperPositionMode','auto');
eval(['print -dpng /Users/Cecilia/TRED46/Satellitbilder/para_interp/ov_',num2str(c),'-',num2str(s),'.png']);
close(gcf)
end

%% tred_surr_sat
surrounding_satellites=tred_surr_sat(532,7,2);
tred_plot_cores(surrounding_satellites{1},surrounding_satellites{2},surrounding_satellites{3})
%% tred_list_cores
%unix('find /Users/Cecilia/TRED46/VirtualSatellite/* >satellites.txt');
available_cores=tred_list_cores;

%% tred_download
[to_download wanted existing]=tred_download(1400:1650);

%% tred_interpolate
interp_sat=tred_interpolate_b(532,8,1700,2);
ind=find(interp_sat{4}<0.5);
%tred_plot_cores(interp_sat{1}(ind,:),interp_sat{2}(ind,:),interp_sat{3}(ind,:))

