% get overall time interval
cd /Users/Cecilia/Data/20090916
tst=irf_time([2009 09 16 07 25 00]);dt=1.5*60;
tint=[tst tst+dt];
%% load fgm and staff data
sc=1:4;
load mBSCR
c_eval('dB?staff=wBSC?;',sc);
c_eval('diB?staff=c_coord_trans(''SR2'',''ISR2'',dB?staff,''cl_id'',?);',sc);
c_eval('diB?fgm=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sc);
c_eval('diB?=c_fgm_staff_combine(diB?fgm,diB?staff);',sc);
%% load efw data
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sc);
c_eval('diE?ib=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sc);
c_eval('diE?ib=irf_edb(diE?ib,diB?);')
%% substitute diE1 with diE1ib since only ib exist where waves are
diE1=diE1ib;
%% get time intervals
sc=1:4;
h=irf_plot('diE?');
c_eval('[t?,~]=get_time(2);',sc); % no bm for C1, only ib
%% refine time intervals
h=irf_plot('diE?');
c_eval('irf_zoom(h,''x'',t?''); [t?,~]=get_time(2);',sc)
%% take out e and b
sc=1:4;
c_eval('Bmatch?=irf_tlim(diB?,t?);',sc)
c_eval('Ematch?=irf_tlim(diE?,t?);',sc)
%% Direction
%for f_highpass
angles=1:3:360;
f_highpass=1:15;
c_eval('[x? y? z? corr_dir? intEdt? Bz? B0? dEk? dEn? Ek? En?]=irf_match_phibe_dir(Bmatch?,Ematch?,angles,f_highpass);',sc);
c_eval('i_dir?=find(corr_dir?(:,1)==max(corr_dir?(:,1)));',sc)
c_eval('direction?=x?(i_dir?,:);',sc)

%% Velocity and density
f_spec=5;
c_eval('[~,i_f?]=min(abs(f_highpass-f_spec));',sc)

n=linspace(0.05,2,50);
v=linspace(00,2000,50);
c_eval('[corr_v?,phi_E?,phi_B?]=irf_match_phibe_v(B0?,Bz?(:,:,i_f?),intEdt?(:,[1 1+i_dir?],i_f?),n,v);',sc)
n_spec=0.25;
c_eval('[~,i_n?]=min(abs(n-n_spec));',sc)
c_eval('i_v?=find(corr_v?(i_n?,:)==min(corr_v?(i_n?,:)));',sc)
c_eval('velocity?=v(i_v?);',sc)
%% Compare correlations for one frequency
plot(angles,corr_dir1,angles,corr_dir2,angles,corr_dir3,angles,corr_dir4)
legend('1','2','3','4')
%% Compare correlations for several frequency
for k=1:4
    h(k)=subplot(2,2,k);
    eval(['pcolor(h(k),corr_dir',num2str(k),')'])
    shading flat
end
linkaxes(h)
colorbar(h)
%% Figure: specify path
thepath='/Users/Cecilia/Research/DF/Pics/3/';

%% Figure: Direction
c_eval('gif_stuff_dir? = irf_match_phibe_vis(''direction'',x?,y?,z?,corr_dir?,intEdt?,Bz?,Ek?,En?);',sc);
c_eval(['imwrite(gif_stuff_dir?.im,gif_stuff_dir?.map,[thepath,cn_time(t?(1)),''-'',cn_time(t?(2)),''_c?_dir.gif''],''DelayTime'',0.01,''LoopCount'',inf);'],sc);
%imwrite(gif_stuff_dir1.im,gif_stuff_dir1.map,[thepath,cn_time(t1(1)),'_c1_dir.gif'],'DelayTime',0.01,'LoopCount',inf);
%% Figure: Velocity
c_eval('gif_stuff_v? = irf_match_phibe_vis(''velocity'',phi_E?,phi_B?(:,[1 i_n?]),v,n(i_n?));',sc);
c_eval('imwrite(gif_stuff_v?.im,gif_stuff_v?.map,[thepath,cn_time(t?(1)),''-'',cn_time(t?(2)),''_c?_v.gif''],''DelayTime'',0.01,''LoopCount'',inf);',sc);

%% Figure: Correlation, density vs velocity
figure; 
for k=sc
hca=subplot(2,2,k); 
c_eval(['axis_handle1=irf_match_phibe_vis(''velocity/density'',hca,n,v,corr_v',num2str(k),');'],k);
c_eval('t1str?=datestr(epoch2date(t?(1)),''dd-mmm-yyyy  HH:MM:SS.FFF'');');
c_eval('t2str?=datestr(epoch2date(t?(2)),''HH:MM:SS.FFF'');');
c_eval('title(hca,[''Correlation C?: '',t1str,'' - '',t2str]);',k)
end
set(gcf,'paperpositionmode','auto','position',[560 260 941 688]);
eval(['print -dpng ',thepath,'vn.png'])

%% Figuer: Compare velocity vectors
scale=1;%[velocity1; velocity2; velocity3; velocity4];
%scale=1;
%scale=[];
cn_plot3d([0 0 0],[direction1; direction2; direction3; direction4],scale,{'C1','C2','C3','C4'},{})
axis square
%% Plot potential
dt=0.035;
h=irf_plot(4);isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,{phi_E3(:,[1 i_v3]),phi_E4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{En3(:,[1 i_v3]),En4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{Ek3(:,[1 i_v3]),Ek4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{irf_tlim(eVTe3,intEdt3(1,1),intEdt3(end,1)),irf_tlim(eVTe4,intEdt(1,1),intEdt(end,1))},'comp')
irf_zoom(h,'x',[intEdt3(1,1),intEdt3(end,1)])