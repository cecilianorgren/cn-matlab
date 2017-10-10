%% Andrej's DF
tint_iso='2009-09-16T07:24:36.300002Z 2009-09-16T07:26:36.300002Z';
tint=irf_time(tint_iso,'iso2epoch')';

%% download from CAA
if 0
prod{1}='C?_CP_EFW_L2_EB';
prod{2}='C?_CP_EFW_L2_E3D_INERT';
prod{3}='C?_CP_FGM_FULL_ISR2';
for k=1:3; caa_download(tint,prod{k}); end
end
%% download staff, done 
if 1
tst=irf_time([2009 9 16 7 10 0]);dt=30*60;
c_get_batch(tst,dt,'vars','bsc','noproc');
end
%% load data
sclist=1:4;
c_eval('diB?fgm=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
c_eval('diE?ib=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sclist);

%% make combined B data from STAFF and FGM
sclist=1:4;
%load mBSCR;
c_eval('diB?staff=c_coord_trans(''SR2'',''ISR2'',wBSC?,''cl_id'',?);',sclist);
c_eval('diB?staffmod=irf_tappl(diB?staff,''*(-1)'');',sclist);
c_eval('diB?=c_fgm_staff_combine(diB?fgm(:,1:4),diB?staff(:,1:4));',sclist);
c_eval('diB?mod=c_fgm_staff_combine(diB?fgm(:,1:4),diB?staffmod(:,1:4));',sclist);
%% make gsm
c_eval('gsmE?=c_coord_trans(''dsc'',''gsm'',diE?,''CL_ID'',?);',sclist);
c_eval('diE?ib=irf_edb(diE?ib,diB?fgm);',sclist);
c_eval('gsmE?ib=c_coord_trans(''isr2'',''gsm'',diE?ib,''CL_ID'',?);',sclist);
c_eval('gsmB?fgm=c_coord_trans(''isr2'',''gsm'',diB?fgm,''CL_ID'',?);',sclist);
c_eval('gsmB?=c_coord_trans(''isr2'',''gsm'',diB?,''CL_ID'',?);',sclist);

%% easy plot
figure;irf_plot('diB?fgm');
figure;irf_plot('diE?');
figure;irf_plot('diE?ib');

%% compare E and B
%h=irf_plot(6);
%isub=0;
title_str='';
print_str='';
fig=figure(99);
set(fig,'position',[560 44 600 900])
setupfigure;
h=irf_plot(9);
isub=1;
if 1 % B_X fgm
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1fgm(:,[1 2]),diB2fgm(:,[1 2]),diB3fgm(:,[1 2]),diB4fgm(:,[1 2])},'comp')
    ylabel(hca,'B_X')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Y fgm
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1fgm(:,[1 3]),diB2fgm(:,[1 3]),diB3fgm(:,[1 3]),diB4fgm(:,[1 3])},'comp')
    ylabel(hca,'B_Y')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Z fgm
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1fgm(:,[1 4]),diB2fgm(:,[1 4]),diB3fgm(:,[1 4]),diB4fgm(:,[1 4])},'comp')
    ylabel(hca,'B_Z')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_X staff
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1staff(:,[1 2]),diB2staff(:,[1 2]),diB3staff(:,[1 2]),diB4staff(:,[1 2])},'comp')
    ylabel(hca,'B_X')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Y staff
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1staff(:,[1 3]),diB2staff(:,[1 3]),diB3staff(:,[1 3]),diB4staff(:,[1 3])},'comp')
    ylabel(hca,'B_Y')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % B_Z staff
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diB1staff(:,[1 4]),diB2staff(:,[1 4]),diB3staff(:,[1 4]),diB4staff(:,[1 4])},'comp')
    ylabel(hca,'B_Z')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % E_X
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diE1(:,[1 2]),diE2(:,[1 2]),diE3(:,[1 2]),diE4(:,[1 2])},'comp')
    ylabel(hca,'E_X')
    irf_legend(hca,{'C1','C2','C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 1 % E_Y
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{diE1(:,[1 3]),diE2(:,[1 3]),diE3(:,[1 3]),diE4(:,[1 3])},'comp')
    ylabel(hca,'E_Y')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
if 0 % E_Z
    hca=h(isub); isub=isub+1;
    irf_plot(hca,{E3(:,[1 4]),E4(:,[1 4])},'comp')
    ylabel('E_Z')
    irf_legend(hca,{'C3','C4'},[0.02 0.95])
    grid(hca,'off')
end
irf_zoom(h,'x',tint)
title(h(1),title_str)
eval(['print -dpng /Users/Cecilia/DF/Pics/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_',print_str,'.gif']);

%% tool
sc=1;
tint_iso_zoom_c4='2009-09-16T07:25:38.000002Z 2009-09-16T07:25:39.000002Z';
tint_iso_zoom_c4='2009-09-16T07:25:37.500002Z 2009-09-16T07:25:38.000002Z';
tint_iso_zoom_c4='2009-09-16T07:25:36.700002Z 2009-09-16T07:25:37.200002Z';
tint=irf_time(tint_iso_zoom_c4,'iso2epoch')';
c_eval('E?=irf_tlim(gsmE?ib,tint);',sc)
c_eval('B?=irf_tlim(gsmB?fgm,tint);',sc)
c_eval('B=B?;',sc);
c_eval('E=E?;',sc);

% Direction
f_filt=1; % Hz
[x y z corr phiE Bz Ek En ufEn ufEk]=tool_direction(irf_tlim(B,tint),irf_tlim(E,tint),300,f_filt);
index_dir=find(corr(:,1)==max(corr(:,1)));
x(index_dir,:)

% Velocity
if size(B,2)==4; B=irf_abs(B); end
B0=mean(irf_tlim(B,tint),1); B0=B0(5);
vint=[100,2000];
title_str=['C',num2str(sc,'%0.f'),', GSM,   f_{filt} = ',num2str(f_filt,'%0.f'),...
    'Hz,  khat = [',num2str(x(index_dir,1),'%0.2f'),' ',...
    num2str(x(index_dir,2),'%0.2f'),' ',num2str(x(index_dir,3),'%0.2f'),']'];
[v correl,A,im,map]=tool_velocity(B0,Bz,phiE(:,[1 index_dir+1]),0.08,vint,50,title_str);
index_v=find(correl(:,1)==min(correl(:,1)));
v(index_v)
tint=[Bz(1,1) Bz(end,1)];
imwrite(im,map,['/Users/Cecilia/DF/Pics/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_C',num2str(sc,'%0.f'),'_f',num2str(sc,'%0.f'),'_velocity.gif'],'DelayTime',0.01,'LoopCount',inf);

%% Illustrate
title_str=['C',num2str(sc,'%0.f'),', GSM,   f_{filt} = ',num2str(f_filt,'%0.f'),'Hz'];
[A,im,map]=tool_gif(x,y,z,corr,[phiE(:,1) phiE(:,2:end)*v(index_v)],ufEk,ufEn,Bz,title_str);
tint=[Bz(1,1) Bz(end,1)];
imwrite(im,map,['/Users/Cecilia/DF/Pics/tool_T',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'-',datestr(fromepoch(tint(2)),'HHMMSSFFF'),'_C',num2str(sc,'%0.f'),'_f',num2str(f_filt,'%0.f'),'_.gif'],'DelayTime',0.01,'LoopCount',inf);
