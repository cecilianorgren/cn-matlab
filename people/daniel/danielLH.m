% The event
myEvent = 1;
if myEvent
    cd /Users/Cecilia/Data/BM/20070831/
    t1 = toepoch([2007 08 31 10 19 05.5]);
    t2 = toepoch([2007 08 31 10 19 05.9]);
    tint = [t1 t2];
else 
    cd /Users/Cecilia/Data/DanielsEvent/
    tint = toepoch([2008 04 22 18 07 32;2008 04 22 18 07 36.9])'; % A
end
% irf_match_phibe_dir.m
% irf_match_phibe_v.m
% irf_match_phibe_vis.m
%% Load and prepare data
%load mBSC.mat % diBSC3,4
if myEvent; 
    load mBS
    c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
    c_eval('diB?fgm=c_coord_trans(''gse'',''isr2'',gseB?fgm,''cl_id'',?);',3:4);
    c_eval('diB?fgm=irf_resamp(diB?fgm,dBS?);',3:4);
    c_eval('diB?=c_fgm_staff_combine(diB?fgm(:,1:4),dBS?);',3:4);
    c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);
else
    load mBSC
    c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
    c_eval('diB?fgm=c_coord_trans(''gse'',''isr2'',gseB?fgm,''cl_id'',?);',3:4);
    c_eval('diB?fgm=irf_resamp(diB?fgm,diBSC?);',3:4);
    c_eval('diB?=c_fgm_staff_combine(diB?fgm(:,1:4),diBSC?);',3:4);
    c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);
end



%% Do the analysis
if myEvent
    tint = toepoch([2008 04 22 18 07 33.65;2008 04 22 18 07 33.85])'; % A
    t1 = toepoch([2007 08 31 10 19 05.4]);
    t2 = toepoch([2007 08 31 10 19 06.1]);
    tint = [t1 t2];
else
    tint = toepoch([2008 04 22 18 07 33;2008 04 22 18 07 34])'; % A
    %tint = toepoch([2008 04 22 18 07 35;2008 04 22 18 07 36])'; % B
end


sc = 3;
if myEvent
    c_eval('E = irf_tlim(diE?,tint);',sc);
    c_eval('B = irf_tlim(diB?,tint);',sc);
else
    c_eval('E = irf_tlim(diE?,tint);',sc);
    c_eval('B = irf_tlim(diB?,tint);',sc);
end

% Direction
angles=1:3:360;
f_highpass=20;
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Velocity and density
n=linspace(0.5,10,100); % cc
v=linspace(1,250,100); % km/s
[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
n_pick=4; % cc
n_diff=abs(n-n_pick);
i_n=find(n_diff==min(n_diff));
i_v=find(corr_v(i_n,:)==min(corr_v(i_n,:)));
velocity=v(i_v);

% Figures
if 1;%plotDirection
    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En); 
    dStrPrint = [irf_time(tint(1),'yyyymmddhhmmss') '-' irf_time(tint(2),'yyyymmddhhmmss') '_f' num2str(f_highpass) 'Hz_sc' num2str(sc) '_d'];
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[dStrPrint '.gif'],'DelayTime',0.01,'LoopCount',inf);
end

if 0; %plotVelocity
    %i_n=50; % if more than one densitiy, choose one by specifying index
    gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
    vStrPrint = [irf_time(tint(1),'yyyymmddhhmmss') '-' irf_time(tint(2),'yyyymmddhhmmss') '_f' num2str(f_highpass) 'Hz_sc' num2str(sc) '_v'];
    imwrite(gif_stuff_v.im,gif_stuff_v.map,[vStrPrint '.gif'],'DelayTime',0.01,'LoopCount',inf);
end

if 0; %plotDensityVelocity
    figure; h(1)=subplot(3,1,1); h(2)=subplot(3,1,2); h(3)=subplot(3,1,3);%axes;
    axis_handle1 = irf_match_phibe_vis('velocity/density',h(1),n,v,corr_v);
    axis_handle2 = irf_match_phibe_vis('velocity/density',h(2),n,v,corr_v);
    axis_handle3 = irf_match_phibe_vis('velocity/density',h(3),n,v,corr_v);
    cStrPrint = [irf_time(tint(1),'yyyymmddhhmmss') '-' irf_time(tint(2),'yyyymmddhhmmss') '_f' num2str(f_highpass) 'Hz_sc' num2str(sc) '_c'];
    %print('-dpng',['corr3' '.png'])
    %cn.colormap(axis_handle1.ax,'default')
end