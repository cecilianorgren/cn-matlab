tint = irf.tint('2015-10-16T10:33:24.50Z',0.9); % magnetosphere-magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:33:44.10Z',0.6); % magnetosphere-magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:33:30.00Z',1); % 
%tint = irf.tint('2015-10-16T10:33:26.90Z',0.2); % 
tint = irf.tint('2015-10-16T10:33:24.00Z',1.5); % 
%tint = irf.tint('2015-10-16T10:33:26.30Z',0.6); % 
%tint = irf.tint('2015-10-16T10:33:26.50Z',1); % 
%tint = irf.tint('2015-10-16T10:33:28.00Z',1); % 
%%
doPlot = 1;
makeFigure = 0;
sc = 3; 
E = dslE3brst.tlim(tint);
acB = dmpaB3scm.tlim(tint);
dcB = dmpaB3.tlim(tint);
c_eval('n_loc = mean(ne_psd?.tlim(tint).data);',sc)
c_eval('Te_loc = mean(Te?perp.tlim(tint).data);',sc)

angles = 0:2:359;
f_highpass = 0.5;

% Propagation direction
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir({dcB,acB},E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Velocity and density
mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3.mean

v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
phiB_scale = B0*1e-18/mu0/e/n; % km/s
v = v_approx*linspace(0.1,2.5,30);

[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
i_v=find(corr_v==min(corr_v));
velocity=v(i_v);

if doPlot % plot corr_dir, corr_v, phiE, phiB
    %% Plot results simply
    nrows = 4;
    hca = subplot(nrows,1,1);
    plot(hca,corr_dir)
    
    hca = subplot(nrows,1,2);
    plot(hca,v,corr_v)
    
    hca = subplot(nrows,1,3);
    tsBz = irf.ts_scalar(irf_time(Bz(:,1),'epoch>epochunix'),Bz(:,2));
    tsIntEdt = irf.ts_scalar(irf_time(intEdt(:,1),'epoch>epochunix'),intEdt(:,1+i_dir));
    irf_plot(hca,{tsBz*phiB_scale,tsIntEdt*velocity},'comp')
    irf_legend(hca,{'B_z*scale','\int Edt*v'},[0.98 0.95])
    
    hca = subplot(nrows,1,4);
    irf_plot(hca,{phi_B,phi_E(:,[1 1+i_v])},'comp')
    irf_legend(hca,{'\phi_B','\phi_E'},[0.98 0.95])
end
if makeFigure 
    %% Figures
    savePath = '/Users/Cecilia/Research/2015Oct16/LH/';
    saveName = [irf_ssub('direction_mms?',ic,f_highpass) '_' irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm') '_' irf_time(tint(2),'epochtt>utc_yyyymmddTHHMMSS.mmm') '_fhp!'];
    
    [out,mva_l,mva_v]=irf_minvar(dslE3brst);
    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath saveName '.gif'],'DelayTime',0.01,'LoopCount',inf);
%%
    i_n=1; % if more than one densitiy, choose one by specifying index
    gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n_loc(i_n));
    imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);

    figure; h=axes;
    axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
end