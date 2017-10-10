%tt = get_time(2); tint = irf_time(tt,'epoch>epochtt');
% Overview plot

h = irf_plot(7);

hca = irf_panel('B fg');
c_eval('irf_plot(hca,{dmpaB?.x,dmpaB?.y,dmpaB?.z,dmpaB?.abs},''comp'');',ic)
hca.YLabel.String = {'fg B_{fg,DMPA}','[nT]'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);


%%
flim = 0.1;
angles=1:3:360;
%f_highpass = flim*flh_loc;
f_highpass = 10;
ic = 1;
%csysB
%csycE

c_eval(['E = dslE?brst.tlim(tint);'],ic)
c_eval(['B = B?sc.tlim(tint);'],ic);

% irf_plot({E,B})
% Get maximum/minimum variance direction for the eletric field
[~,mva_l,mva_v]=irf_minvar(E.tlim(tint));


% Propagation direction
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

%% Vizualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
%gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
savePath = '/Users/Cecilia/Research/2015Oct01/';
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'C' num2str(ic) '_' irf_time(tint(1),'epochtt>utc_HH:MM:SS.mmm') '_' num2str(f_highpass,'%.1f') '_dir_aaaaaaa.gif'],'DelayTime',0.01,'LoopCount',inf);
%imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);