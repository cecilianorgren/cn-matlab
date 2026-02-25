% For Alexandros Chasapis
cd /Users/Cecilia/Data/Cluster/20020327
savePath = '/Users/Cecilia/Research/LH2/20020327/';
tint_waves = toepoch([2002 03 27 10 14 02.9;2002 03 27 10 14 03.2])'; 
tint_dl = toepoch([2002 03 27 10 10 00;2002 03 27 10 20 00])'; 

doDownload = 0;
doLoad = 0;
doLookAtData = 0;
doDeriveParameters = 0;
if doDownload
    %%  
    cd /Users/Cecilia/Data/Cluster/20020310/
    scList = 3:4;
    c_eval([,...
    'caa_download(tint_dl,''C?_CP_EFW_L2_E'');',...
    'caa_download(tint_dl,''C?_CP_EFW_L2_E3D_GSE'');',...
    'caa_download(tint_dl,''C?_CP_STA_CWF_GSE'');',...
    'caa_download(tint_dl,''C?_CP_FGM_FULL'');',...
    'caa_download(tint_dl,''C?_CP_PEA_MOMENTS'');',...
    ],scList);
end
if doLoad
    c_eval('gseE? = c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''mat''); gseE?(:,1)=gseE?(:,1)-1/450;',sclist);
    c_eval('gseB?fgm = c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
    c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
    c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?sta,''cl_id'',?);',sclist)
    c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
    c_eval('parTe?MK=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
end
if doDeriveParameters
    units = irf_units;
    c_eval('absB? = irf_abs(gseB?fgm); absB?=absB?(:,[1 5]);',sclist)
    c_eval('flh? = [absB?(:,1) absB?(:,2)*1e-9*units.e/sqrt(units.me*units.mp)/2/pi];',sclist)
    c_eval('absE? = irf_abs(gseE?); wavE? = irf_wavelet(absE?(:,[1 5]));',sclist);    
    c_eval('peaNe?hf=irf_resamp(peaNe?,gseE?);',sclist);
    c_eval('parTe?=[parTe?MK(:,1) irf_convert(parTe?MK(:,2),''MK2eV'')];',sclist);
    c_eval('vte? = [parTe?(:,1) sqrt(2*units.e*parTe?(:,2)/units.me)*1e-3];',sclist) % km/s
    c_eval('fce? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.me/2/pi];',sclist)
    c_eval('re? = irf_multiply(1/2/pi/sqrt(2),vte?,1,fce?,-1);',sclist) % km
end
if doLookAtData
    %%
    h = irf_plot(5,'newfigure');
    hca = irf_panel('Bx');
    irf_plot(hca,{gseB3fgm(:,[1 2]),gseB3(:,[1 2])},'comp')
    ylabel(hca,'B_x [nT]')
    
    hca = irf_panel('By'); 
    irf_plot(hca,{gseB3fgm(:,[1 3]),gseB3(:,[1 3])},'comp')
    ylabel(hca,'B_y [nT]')
    
    hca = irf_panel('Bz'); 
    irf_plot(hca,{gseB3fgm(:,[1 4]),gseB3(:,[1 4])},'comp')
    ylabel(hca,'B_z [nT]')
    
    hca = irf_panel('staB');
    irf_plot(hca,gseB3sta)
    ylabel(hca,'B [nT] (STAFF)')
    irf_legend(hca,{'x','y','z'},[0.98 0.95])
    
    hca = irf_panel('E');
    irf_plot(hca,gseE3)
    ylabel(hca,'E [mV/m]')
    irf_legend(hca,{'x','y','z'},[0.98 0.95])
    
    if 0
    hca = irf_panel('E wavelet');
    irf_spectrogram(hca,wavE3);
    hold(hca,'on');
    irf_plot(hca,flh3,'color','white')
    ylabel(hca,'E')
    end
    
    irf_zoom(h,'x',tint_waves)
    irf_zoom(h,'y')
end
if 0
    %%
    h = irf_plot(2,'newfigure');
    
    hca = irf_panel('B');
    irf_plot(hca,gseB3)
    ylabel(hca,'B [nT] (GSE)')
    irf_legend(hca,{'x','y','z'},[0.98 0.95])
    
    hca = irf_panel('E');
    irf_plot(hca,gseE3)
    ylabel(hca,'E [mV/m] (GSE)')
    irf_legend(hca,{'x','y','z'},[0.98 0.95])
    
    % add lower hybrid interval
    tint_flh = tint_waves(1)-0.18;
    c_eval('flh_loc = irf_resamp(flh?,tint_flh,''nearest'');',sc); flh_loc = flh_loc(2);
    tint_flh = tint_flh + [-1 1]/flh_loc;
    irf_pl_mark(h,tint_flh,[1 0.2 0]);

    % add lower hybrid interval
    tint_flh = tint_waves(1)-0.4;
    c_eval('flh_loc = irf_resamp(flh?,tint_flh,''nearest'');',sc); flh_loc = flh_loc(2);
    tint_flh = tint_flh + [-1 1]/flh_loc;
    irf_pl_mark(h,tint_flh,[1 0.7 0]);
    
    % add lower hybrid interval
    tint_flh = tint_waves(1)-0.8;
    c_eval('flh_loc = irf_resamp(flh?,tint_flh,''nearest'');',sc); flh_loc = flh_loc(2);
    tint_flh = tint_flh + [-1 1]/flh_loc;
    irf_pl_mark(h,tint_flh,[1 0.7 0]);
    
    % add lower hybrid interval
    tint_flh = tint_waves(1);
    c_eval('flh_loc = irf_resamp(flh?,tint_flh,''nearest'');',sc); flh_loc = flh_loc(2);
    tint_flh = tint_flh + [-1 1]/flh_loc;
    irf_pl_mark(h,tint_flh,[1 0.7 0]);
    
    
    irf_pl_mark(h,tint_waves,[0.8 0.4 0.8]);
    
    legend(hca,'crossing','1/<f_{LH}>','1/<f_{LH}>','1/<f_{LH}>','1/<f_{LH}>','location','southeast')
    
    
    %irf_zoom(h,'x',tint_dl+[3.75*60 -5.75*60])
    irf_zoom(h,'x',tint_waves+[-2 1])
    irf_zoom(h,'y')
end

sc = 3;
if 0
    tint = tint_waves;
elseif 0
    %c_eval('irf_plot({irf_tlim(gseE?,tint_waves),irf_tlim(gseB?,tint_waves)})',sc)
    tint = get_time(2); tint = torow(tint);
else
    tint = toepoch([2002 03 27 10 14 02.85;2002 03 27 10 14 03])'; 
    tint = toepoch([2002 03 27 10 14 02.65;2002 03 27 10 14 02.8])'; 
end

csys = 'gse';
flim = 0.7;
v=linspace(20,800,30);

% Obtain correlation parameters
% Assume data is loaded.
% Also assume the spacecraft is given and that the timeinterval is defined
% in som way beforehand. 
c_eval('re_loc = irf_resamp(re?,tint(1),''nearest'');',sc); re_loc = re_loc(2);
c_eval('flh_loc = irf_resamp(flh?,tint(1),''nearest'');',sc); flh_loc = flh_loc(2);
angles=1:3:360;
f_highpass = flim*flh_loc;
c_eval(['E = irf_tlim(' csys 'E?,tint);'],sc)
c_eval(['B = irf_tlim(' csys 'B?,tint);'],sc);

% Propagation direction
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Find error
Cerror = numel(find(corr_dir>max(corr_dir)*0.95))*3/2;

% Get maximum/minimum variance direction for the eletric field
[~,mva_l,mva_v]=irf_minvar(irf_tlim(E,tint));

% Velocity and density
c_eval('n_loc = irf_resamp(peaNe?hf,tint(1),''nearest'');',sc); n_loc = n_loc(2);

% Approximate v range
mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3
v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
v = v_approx*linspace(0.5,1.5,30);

[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
i_v=find(corr_v(1,:)==min(corr_v(1,:)));
velocity=v(i_v);

irf_match_phibe_vis('best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc,max(phi_B(:,2))/max(Bz(:,2)));
fig = gcf;
hca = fig.UserData.subplot_handles;
irf_legend(hca,{'',['B_0=' num2str(B0,'%.1f') ' nT  n=' num2str(n_loc,'%.1f') ' cm^{-3}']},[0.99 0.1])
% cn.print([ 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' irf_time(tint(2),'epoch>utc') '_fhighpass' num2str(flim,'%.1f') 'flh_best']);
%% Visualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'C' num2str(sc) '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_dir_.gif'],'DelayTime',0.01,'LoopCount',inf);
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath 'C' num2str(sc) '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);


%% Add the event and save data to TimeTable
comment = 'very good!';
toSave = 'single';
tool.add_event
end