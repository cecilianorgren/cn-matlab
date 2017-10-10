%% Choose spacecraft
sc=3:4;

t1=[2007 09 02 15 47 29.30]; t2=[2007 09 02 15 47 30.80]; tint=toepoch([t1;t2]);
%% Fields to match
sc=4;
c_eval('Bmatch=irf_tlim(diB?staff,t);',sc)
c_eval('Ematch=irf_tlim(diE?,t);',sc)
%% Direction
angles=1:3:360;
f_highpass=6;
c_eval('[x? y? z? corr_dir? intEdt? Bz? B0? dEk? dEn? Ek? En?]=irf_match_phibe_dir(Bmatch,Ematch,angles,f_highpass);',sc);
c_eval('i_dir?=find(corr_dir?(:,1)==max(corr_dir?(:,1)));',sc)
c_eval('direction?=x?(i_dir?,:);',sc)

% Velocity and density
n=linspace(0.05,1,500);
v=linspace(00,1000,500);
c_eval('[corr_v?,phi_E?,phi_B?]=irf_match_phibe_v(B0?,Bz?,intEdt?(:,[1 1+i_dir?]),n,v);',sc)
n_spec=0.5;
c_eval('[~,i_n?]=min(abs(n-n_spec));',sc)
c_eval('i_v?=find(corr_v?(i_n?,:)==min(corr_v?(i_n?,:)));',sc)
c_eval('velocity?=v(i_v?);',sc)

%% Figures
c_eval('x=x?;y=y?;z=z?;corr_dir=corr_dir?;intEdt=intEdt?;Bz=Bz?;Ek=Ek?;En=En?;phi_E=phi_E?;phi_B=phi_B?;',4)
thepath='/Users/Cecilia/Research/DF/Pics/';
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En);      
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[thepath,'mygif_dir4.gif'],'DelayTime',0.01,'LoopCount',inf);

%i_n=50; % if more than one densitiy, choose one by specifying index
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E4,phi_B4(:,[1 i_n4]),v,n(i_n4));
imwrite(gif_stuff_v.im,gif_stuff_v.map,[thepath,'mygif_v3.gif'],'DelayTime',0.01,'LoopCount',inf);

figure; %h1=subplot(1,2,1);axis_handle1 = irf_match_phibe_vis('velocity/density',h1,n,v,corr_v3);
h2=subplot(1,2,2);axis_handle2 = irf_match_phibe_vis('velocity/density',h2,n,v,corr_v4);
%h3=subplot(1,3,3);axis_handle3 = irf_match_phibe_vis('velocity/density',h3,n,v,corr_v4-corr_v3);
%% Plot potential
dt=0.035;
h=irf_plot(4);isub=1;
hca=h(isub);isub=isub+1;irf_plot(hca,{phi_E3(:,[1 i_v3]),phi_E4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{En3(:,[1 i_v3]),En4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{Ek3(:,[1 i_v3]),Ek4(:,[1 i_v4])},'comp','dt',[0 dt])
hca=h(isub);isub=isub+1;irf_plot(hca,{irf_tlim(eVTe3,intEdt3(1,1),intEdt3(end,1)),irf_tlim(eVTe4,intEdt(1,1),intEdt(end,1))},'comp')
irf_zoom(h,'x',[intEdt3(1,1),intEdt3(end,1)])