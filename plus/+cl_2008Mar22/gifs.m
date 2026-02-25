em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
tint_out = toepoch([2008 04 22 17 35 30;2008 04 22 17 37 30])';
tint = em_tint;
tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 34])'+7.5+[0.3 -0.3];
flim = 1;
    tool.single_event   
%% Figures
    [out,mva_l,mva_v]=irf_minvar(gsmE3);
    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,['gif_dir_dr2' num2str(flim) '.gif'],'DelayTime',0.01,'LoopCount',inf);
%%
    i_n=1; % if more than one densitiy, choose one by specifying index
    gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n_loc(i_n));
    imwrite(gif_stuff_v.im,gif_stuff_v.map,'gif_v_dr1.gif','DelayTime',0.01,'LoopCount',inf);

    figure; h=axes;
    axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);

%% Figures
    [out,mva_l,mva_v]=irf_minvar(gsmE3);
    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);

    i_n=1; % if more than one densitiy, choose one by specifying index
    gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n_loc(i_n));
    imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);

    figure; h=axes;
    axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
