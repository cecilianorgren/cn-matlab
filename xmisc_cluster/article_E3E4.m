diB3=c_coord_trans('gsm','dsi',gsmB3,'cl_id',3);
diB4=c_coord_trans('gsm','dsi',gsmB4,'cl_id',4);
%%
[diE3_20,ang3_20]=irf_edb(irf_tlim(diE3,tint(1),tint(2)),diB3,20);
[diE4_20,ang4_20]=irf_edb(irf_tlim(diE4,tint(1),tint(2)),diB4,20);
[diE3_15,ang3]=irf_edb(irf_tlim(diE3,tint(1),tint(2)),diB3,15);
[diE4_15,ang4]=irf_edb(irf_tlim(diE4,tint(1),tint(2)),diB4,15);
%%
gsmE3_15=c_coord_trans('dsi','gsm',diE3_15,'cl_id',3);
gsmE4_15=c_coord_trans('dsi','gsm',diE4_15,'cl_id',4);
gsmE3_20=c_coord_trans('dsi','gsm',diE3_20,'cl_id',3);
gsmE4_20=c_coord_trans('dsi','gsm',diE4_20,'cl_id',4);
gsmE3_25=c_coord_trans('dsi','gsm',diE3_25,'cl_id',3);
gsmE4_25=c_coord_trans('dsi','gsm',diE4_25,'cl_id',4);
%%
gseE3_15=c_coord_trans('dsi','gse',diE3_15,'cl_id',3);
gselE4_15=c_coord_trans('dsi','gse',diE4_15,'cl_id',4);
%%
figure;h=irf_plot(6);
irf_plot(h(1),diE3_20(1:fix(end/2),[1 4]));
irf_plot(h(2),diE3_25(1:fix(end/2),[1 4]));
irf_plot(h(3),[diE3_20(1:fix(end/2),1) tocolumn(ang3(1:fix(end/2)))]);
irf_plot(h(4),diE4_20(1:fix(end/2),[1 4]));
irf_plot(h(5),diE4_25(1:fix(end/2),[1 4]));
irf_plot(h(6),[diE4_25(1:fix(end/2),1) tocolumn(ang4(1:fix(end/2)))]);
irf_zoom(h,'x',[toepoch([2007 08 31 10 13 00]) toepoch([2007 08 31 10 24 00])])
%%
figure;h=irf_plot(6);
irf_plot(h(1),gsmE3_20(1:fix(end/2),[1 3]));
irf_plot(h(2),gsmE3_25(1:fix(end/2),[1 3]));
irf_plot(h(3),[gsmE3_20(1:fix(end/2),1) tocolumn(ang3(1:fix(end/2)))]);
irf_plot(h(4),gsmE4_20(1:fix(end/2),[1 3]));
irf_plot(h(5),gsmE4_25(1:fix(end/2),[1 3]));
irf_plot(h(6),[gsmE4_25(1:fix(end/2),1) tocolumn(ang4(1:fix(end/2)))]);
irf_zoom(h,'x',[toepoch([2007 08 31 10 13 00]) toepoch([2007 08 31 10 24 00])])
%%
difromgseE4=c_coord_trans('gse','dsi',gseE4,'cl_id',4);
difromgseE3=c_coord_trans('gse','dsi',gseE3,'cl_id',3);
%%
for k=1:8    
    angles=k*5;
    n_nan(k)=max(size(find(isnan(irf_edb(irf_tlim(diE4,tint(1),tint(2)),irf_tlim(diB4,tint(1),tint(2)),angles))==1)));
end
%%
the_angle=18;
[a,b]=irf_edb(irf_tlim(diE4,tint(1),tint(2)),diB3,the_angle);
figure;h=irf_plot({irf_tlim(a(:,[1 4]),tint(1),tint(2)),[a(:,1) ang3 ang3*0+the_angle]});
title(h(1),['Limiting angle = ',num2str(the_angle)])
set(gcf,'PaperPositionMode','auto');
eval(['print -dpdf constructing_limitingangle',num2str(the_angle),'.pdf']);