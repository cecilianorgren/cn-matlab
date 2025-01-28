tint = irf.tint('2015-10-30T05:15:44.00Z/2015-10-30T05:15:48.00Z');
tint = irf.tint('2015-10-30T05:15:46.10Z/2015-10-30T05:15:46.80Z');
tintOverview = irf.tint('2015-10-30T05:15:44.00Z/2015-10-30T05:15:50.00Z');
tint = irf.tint('2015-10-30T05:15:46.60Z/2015-10-30T05:15:46.70Z');
tint = irf.tint('2015-10-30T05:15:46.50Z/2015-10-30T05:15:46.80Z'); % mms2
tint = irf.tint('2015-10-30T05:15:46.55Z/2015-10-30T05:15:46.65Z'); % mms2
tint = irf.tint('2015-10-30T05:15:45.70Z/2015-10-30T05:15:46.00Z'); % mms2
tint = irf.tint('2015-10-30T05:15:46.70Z/2015-10-30T05:15:46.90Z'); % mms2
tint = irf.tint('2015-10-30T05:15:46.70Z/2015-10-30T05:15:47.10Z'); % mms1
%%
doPlot = 1;
makeFigure = 0;
ic = 1; 
c_eval('E = gseE?.tlim(tint);',ic)
c_eval('acB = gseB?scm.tlim(tint);',ic)
c_eval('dcB = gseB?.tlim(tint);',ic)
c_eval('n_loc = mean(ne?.tlim(tint).data);',ic)
c_eval('Te_loc = mean(facTe?.tlim(tint).yy.data + facTe?.tlim(tint).zz.data)*0.5;',ic)

angles = 0:3:359;
f_highpass = 30;5*8192/450;40;

% Propagation direction
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir({dcB,acB},E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Velocity and density
mu0=4*pi*1e-7; e=1.6e-19; n=n_loc*1e6; % density in #/m^3.mean

v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
phiB_scale = B0*1e-18/mu0/e/n; % km/s, This is a fix scale for given time interval
v = v_approx*linspace(0.1,2.5,30);

[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v);
i_v=find(corr_v==min(corr_v));
velocity=v(i_v);

% make independent estimate of phi from B
lowpassB = f_highpass;
highpassB = f_highpass;
c_eval('fsDFG = 1/(gseB?.time(2)-gseB?.time(1));',ic)
c_eval('fsSCM = 1/(gseB?scm.time(2)-gseB?scm.time(1));',ic)
c_eval('BDC = irf_filt(gseB?,0,lowpassB,fsDFG,3);',ic)
c_eval('BAC = irf_filt(gseB?scm,1,0,fsSCM,3);',ic)
 c_eval('locPhiScale = BDC.abs*1e-9/(1e6*e*4*pi*1e-7)/ne?;',ic)        
PhiB = BAC.dot(BDC.resample(BAC))/BDC.abs*1e-9*locPhiScale.resample(BAC);


if doPlot % plot corr_dir, corr_v, phiE, phiB
    %% Plot results
    [h1,h2] = initialize_combined_plot(12,2,1,0.5,'vertical');
    
    % Plot time series
    % Overview
    hca = irf_panel('DC B');
    irf_plot(hca,{BDC.x,BDC.y,BDC.z,BDC.abs},'comp')
    hca.YLabel.String = 'B_{DC}';
    irf_legend(hca,{['f <  ' num2str(lowpassB) ' Hz']},[0.02 0.95],'color',[0 0 0])
    
    hca = irf_panel('AC B');
    irf_plot(hca,{BAC.x,BAC.y,BAC.z,BAC.abs},'comp')
    hca.YLabel.String = 'B_{AC}';
    irf_legend(hca,{['f >  ' num2str(highpassB) ' Hz']},[0.02 0.95],'color',[0 0 0])
    
    hca = irf_panel('n');
    c_eval('irf_plot(hca,ne?);',ic)
    hca.YLabel.String = 'n_e';
    
    hca = irf_panel('PhiB');
    c_eval('irf_plot(hca,{PhiB},''comp'');',ic); irf_legend(hca,{'\phi_B'},[0.02 0.95])
    %c_eval('irf_plot(hca,{irf.,PhiB},''comp'');',ic); irf_legend(hca,{'|\phi_B|','\phi_B'},[0.02 0.95])
    hca.YLabel.String = '\phi_B';
    
    
    hca = irf_panel('T');
    c_eval('irf_plot(hca,{facTe?.xx,0.5*(facTe?.yy+facTe?.zz)},''comp'');',ic)
    hca.YLabel.String = 'T_e';
    irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.95])
    
    hca = irf_panel('f LH');
    c_eval('irf_plot(hca,{flh?},''comp'')',ic)
    hca.YLabel.String = 'f_{LH}';
    
    delete(irf_panel('delete'))
    %hca = irf_panel('DFG B');
    %irf_plot(hca,{dcB.x,dcB.y,dcB.z,dcB.abs},'comp')
    %hca.YLabel.String = 'B_{DFG}';
    
    hca = irf_panel('wave Bz');
    irf_plot(hca,Bz,'comp')
    hca.YLabel.String = 'B_{SCM}';
    
    hca = irf_panel('E');
    irf_plot(hca,E)
    hca.YLabel.String = 'E';    
    
    hca = irf_panel('Ek');
    irf_plot(hca,{Ek(:,[1 1+i_dir]),dEk(:,[1 1+i_dir])},'comp')
    hca.YLabel.String = 'E_k (mV/m)';
    irf_legend(hca,{'E_k','\delta E_k'},[0.98 0.95])
    irf_legend(hca,{['k = [' num2str(direction,'%.2f') '] (GSE)']},[0.02 0.95],'color',[0 0 0])
    
    
    hca = irf_panel('Phi 0');
    tsBz = irf.ts_scalar(irf_time(Bz(:,1),'epoch>epochTT'),Bz(:,2));
    tsIntEdt = irf.ts_scalar(irf_time(intEdt(:,1),'epoch>epochTT'),intEdt(:,1+i_dir));
    irf_plot(hca,{tsBz*phiB_scale,tsIntEdt*velocity},'comp')
    hca.YLabel.String = '\phi (V)';
    irf_legend(hca,{'B_z*B_0/en\mu_0','\int Edt*v'},[0.98 0.95])
    irf_legend(hca,{['B_0 = ' num2str(B0,'%.1f') ' nT, n = ' num2str(n_loc,'%.1f') ' cc']},[0.02 0.95],'color',[0 0 0])
    
    hca = irf_panel('Phi');
    irf_plot(hca,{phi_B,phi_E(:,[1 1+i_v])},'comp')
    hca.YLabel.String = '\phi (V)';
    irf_legend(hca,{'\phi_B','\phi_E'},[0.98 0.95])
    
    h1(1).Title.String = irf_ssub('MMS ?',ic);
    irf_zoom(h1(1:6),'x',tintOverview)
    irf_zoom(h1(8:end),'x',tint)
    
    irf_zoom(h1(1:6),'y')    
    irf_zoom(h1(8:end),'y')    
    irf_plot_zoomin_lines_between_panels(h1(6),h1(8))
    add_length_on_top(h1(8),velocity,1)
    
    % Plot correlation results
    isub = 1;
    hca = h2(isub); isub = isub + 1;
    plot(hca,angles,corr_dir,angles(i_dir),corr_dir(i_dir),'ro')
    axes(hca)
    text(angles(i_dir),0,['k = [' num2str(direction,'%.2f') '] (GSE)'],'horizontalalignment','center','color','r')
    hca.YLabel.String = 'C_{direction}';
    hca.XLabel.String = 'Angle (deg)';
    
    hca = h2(isub); isub = isub + 1;
    plot(hca,v,corr_v)
    hca.YLabel.String = 'C_{velocity}';
    hca.XLabel.String = 'Velocity (km/s)';
    
    
end
if makeFigure 
    %% Figures
    savePath = '/Users/Cecilia/Research/2015Oct30/LH/';
    saveName = [irf_ssub('direction_mms?',ic,f_highpass) '_' irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm') '_' irf_time(tint(2),'epochtt>utc_yyyymmddTHHMMSS.mmm') '_fhp!'];
    
    [out,mva_l,mva_v]=irf_minvar(gseE3);
    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath saveName '.gif'],'DelayTime',0.01,'LoopCount',inf);
%%
    i_n=1; % if more than one densitiy, choose one by specifying index
    gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n_loc(i_n));
    imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);

    figure; h=axes;
    axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
end