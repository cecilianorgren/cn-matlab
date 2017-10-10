cd /Users/Cecilia/Data/Cluster/20070417/
savePath = '/Users/Cecilia/Research/LH2/20070417/';
sclist=1:4;

%load tt20070417

% Electric field
c_eval('diE? = c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
c_eval('diE?(:,1) = diE?(:,1)-1/450;') % Adjust for time error in CAA data.

% Magnetic field
c_eval('diB?fgm = c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);
c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
c_eval('diB?sta=c_coord_trans(''GSE'',''ISR2'',gseB?sta,''cl_id'',?);',sclist);
c_eval('diB?=c_fgm_staff_combine(diB?fgm(:,1:4),diB?sta,''cl_id'',?);',sclist)
c_eval('absB? = irf_abs(diB?fgm); absB? = absB?(:,[1 5]);',sclist)    

% Electron moments
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?MK=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?=[parTe?MK(:,1) irf_convert(parTe?MK(:,2),''MK2eV'')];',sclist);
%c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('peaNe?hf = irf_resamp(peaNe?,diE?);',sclist);
c_eval('e_deflux?=c_caa_var_get(''Data__C?_CP_PEA_PITCH_SPIN_DEFlux'',''mat'');',sclist);

% Ion moments
c_eval('gseVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',[1 3:4]);
c_eval('gsmVi? = c_coord_trans(''ISR2'',''GSM'',gseVi?,''cl_id'',?);',[1 3]);


% Position
c_eval('diR? = c_caa_var_get(''sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);

% Make 3D electric field data
c_eval('[diE?,ang?]=irf_edb(diE?,diB?fgm);',sclist);

% Transform fields to GSM coordinates
c_eval('gsmE? = c_coord_trans(''ISR2'',''GSM'',diE?,''cl_id'',?);',sclist);
c_eval('gsmB? = c_coord_trans(''ISR2'',''GSM'',diB?,''cl_id'',?);',sclist);


% Derived quantities
units = irf_units;

% Electron thermal velocity
c_eval('vte? = [parTe?(:,1) sqrt(2*units.e*parTe?(:,2)/units.me)*1e-3];',sclist) % km/s

% Electron cyclotron frequency
c_eval('fce? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.me/2/pi];',sclist)
    
% Lower hybrid frequency
c_eval('flh? = [absB?(:,1) absB?(:,2)*1e-9*units.e/sqrt(units.me*units.mp)/2/pi];',sclist)
    
% Electron gyroradius
c_eval('re? = irf_multiply(1/2/pi/sqrt(2),vte?,1,fce?,-1);',sclist) % km

% Integrate E
c_eval('intdiE? = irf_integrate(diE?);',sclist) % mVs/m

% Make power spectrograms
nfft = 1024;
sfreq = 450;
overlap = 0.50;
c_eval('fftE? = irf_powerfft(diE?,nfft,sfreq,[overlap]);',sclist)
c_eval('fftB? = irf_powerfft(diB?,nfft,sfreq,[overlap]);',sclist)

%% Overview plot for single SC 
sc=2;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(7);
isub = 1;

if 1 % Electron differential energy flux
    hca = h(isub); isub=isub+1;
    c_eval('deflux = e_deflux?;',sc)
    
    specrec.t = deflux.t;
    specrec.f = deflux.dep_x{2}.data(1,:)'*1e3; % eV
    specrec.p = {squeeze(nansum(deflux.data,2))};
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f)
    hca.YScale = 'log';
    hca.YTick = [0.01 0.1 1 10 100 1000 10000];
    hca.YLabel.String = 'E_e [eV]';
    hc = colorbar('peer',hca);
    hc.Label.String = 'keV/cm^2 s sr keV';
    
    %irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dEF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    %hca.YLabel.String = 'E [nT]';
end
if 1 % n (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('peaNe?',sc));
    hca.YLabel.String = 'n_e [cc]';    
end
if 1 % Te (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('parTe?',sc));
    hca.YLabel.String = 'T_{e,||} [eV]';    
end

if 0 % B ISR2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 0 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(nT)^2/Hz';
    hca.NextPlot = 'add';
    c_eval('irf_plot(hca,flh?)',sc);
    c_eval('irf_plot(hca,fce?)',sc);
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';
end
if 1 % B power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftB?.t,fftB?.p{1},fftB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(mV/m)^2/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end
if 0
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftEB?.t,fftEB?.p{1},fftEB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = 'mV/m/nT/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end


title(h(1),irf_ssub('C?',sc))
irf_zoom(h,'x',tint)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])
%% Overview plot for 4 SC, fields
sc=1:4;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(8);
isub = 1;

if 1 % B1
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diB1);
    hca.YLabel.String = 'B1 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diB2);
    hca.YLabel.String = 'B2 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diB3);
    hca.YLabel.String = 'B3 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diB4);
    hca.YLabel.String = 'B4 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E1
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diE1);
    hca.YLabel.String = 'E1 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diE2);
    hca.YLabel.String = 'E2 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diE3);
    hca.YLabel.String = 'E3 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,diE4);
    hca.YLabel.String = 'E4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end


title(h(1),'')
irf_zoom(h,'x',tint)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])
%% Overview plot for 2-4 SC, fields (E,V,Vi)
sc=1:4;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(9);
isub = 1;

if 1 % B2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB1);
    hca.YLabel.String = 'B1 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB2);
    hca.YLabel.String = 'B2 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB3);
    hca.YLabel.String = 'B3 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB4);
    hca.YLabel.String = 'B4 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE2);
    hca.YLabel.String = 'E2 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE3);
    hca.YLabel.String = 'E3 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE4);
    hca.YLabel.String = 'E4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % Vi1
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmVi1);
    hca.YLabel.String = 'Vi1 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % Vi3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3);
    hca.YLabel.String = 'Vi4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end


title(h(1),'')
irf_zoom(h,'x',tint)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
%eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])

%% Make gifs of a few time intervals for all spacecraft.
% make timetable
% tt20070417 = irf.TimeTable;
sc = 2;
%irf_plot(irf_ssub('diE?',sc));
%tint = get_time(2); tint = torow(tint);
%tt20070417=add(tt20070417,tint);
tint = toepoch([2007 04 17 15 32 50.70; 2007 04 17 15 32 51.00])';
%%
timeStr = datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff')

angles=1:3:360;
f_highpass=50;
flh_loc=f_highpass;
f_filt_loc=f_highpass;
c_eval('E = irf_tlim(diE?,tint);',sc)
c_eval('B = irf_tlim(diB?,tint);',sc);
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Get maximum variance direction for the eletric field
[~,mva_l,mva_v]=irf_minvar(irf_tlim(E,tint));
%% Velocity and density
n=linspace(0.1,5,100);
v=linspace(20,1000,30);
[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
i_n=34;
i_v=find(corr_v(i_n,:)==min(corr_v(i_n,:)));
velocity=v(i_v);

%% Figures
figure; h=axes;
axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
 
%%       
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,mva_l,mva_v,f_highpass);      
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath irf_time(tint(1),'epoch2iso') '_' num2str(f_highpass) '_dir_.gif'],'DelayTime',0.01,'LoopCount',inf);

%%
i_n=34;
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath irf_time(tint(1),'epoch2iso') '_' num2str(f_highpass) '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);
 
%% Moving average
sc = 2; 
tint = toepoch([2007 04 17 15 28 20; 2007 04 17 15 39 20])';
tint = toepoch([2007 04 17 15 32 00; 2007 04 17 15 36 00])';
tint = toepoch([2007 04 17 15 35 00; 2007 04 17 15 35 15])';
%%
tint = toepoch([2007 04 17 15 34 21.5; 2007 04 17 15 34 23.5])';
tint = toepoch([2007 04 17 15 34 19.0; 2007 04 17 15 34 21.4])';
tint = toepoch([2007 04 17 15 32 00; 2007 04 17 15 36 00])';
%%
ffiltlim = 0.8;
nTlh = 3; 
vrange = [50 800];
c_eval('[mB,cmaxs,mE,mintEdt,ks,ns,Bs,B0s,cs,mva_v1,mva_v2,mva_v3,mva_l123,EB,tints,ffilts,vs] = tool.moving_average(irf_tlim(diE?,tint),irf_tlim(diB?,tint),irf_tlim(flh?,tint),ffiltlim,nTlh,peaNe?hf,vrange);',sc);
%% Plot moving average
h = irf_plot(9);
isub = 1;

if 1 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_tlim(diE2,tint));
    hca.YLabel.String = 'E2 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E maximum variance direction
    hca = h(isub); isub=isub+1;
    irf_plot(hca,mva_v1)    
    hca.YLabel.String = 'Maximum \newline variance \newline direction';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % k from matching
    hca = h(isub); isub=isub+1;
    irf_plot(hca,ks)  
    hca.YLabel.String = 'k';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;
    irf_plot(hca,cmaxs)  
    hca.YLabel.String = 'C_{max}';      
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;   
    specrec.t = tints(:,1);
    specrec.p = cs;
    specrec.f = 1:3:360';
    irf_spectrogram(hca,specrec) 
    hca.YLabel.String = '\theta';
    hca.CLim = [-1 1];
    hc = colorbar('peer',hca);
    hc.YLim = [-1 1];
    hc.YLabel.String = 'C';
    cmap = irf_colormap('poynting');
    colormap(hca,cmap)
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,mE)    
    hca.YLabel.String = '<|\delta E|>';    
end

if 1 % frequencies 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,ffilts); hold(hca,'on');
    c_eval('flhplot = irf_tlim(flh?,tint);',sc)
    irf_plot(hca,flhplot); hold(hca,'off');
    hca.YLabel.String = 'f [Hz]';
    irf_legend(hca,{'f_{filter}','f_{LH}'},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end
if 1 % frequencies 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,mva_l123); 
    hca.YLabel.String = 'l_{MVA,E}';
    irf_legend(hca,{'max','inter','min'},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end
if 1 % velocitiess 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,vs);     
    clim = 0.6; 
    ind_plot = find(cmaxs(:,2)>clim);
    hold(hca,'on')
    irf_plot(hca,vs(ind_plot,:),'r*');
    hca.YLabel.String = 'v [km/s]';
    irf_legend(hca,{'','',['v(C>' num2str(clim) ')']},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end

title(h(1),['nTlh = ' num2str(nTlh) ', f_{filt} = ' num2str(ffiltlim) 'f_{LH}'])
irf_zoom(h,'x',tint)
irf_plot_axis_align

%% Make turning gif of certain time interval
ind = 22;

for ind = [13 17 20 59];22;[7 11 14 16 19 26 29 31];
    tint = tints(ind,2:3);
    angles=1:3:360;
    f_highpass=ffilts(ind,2);
    c_eval('E = irf_tlim(diE?,tint);',sc)
    c_eval('B = irf_tlim(diB?,tint);',sc);
    [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
    i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
    direction=x(i_dir,:);

    % Get maximum variance direction for the eletric field
    [~,mva_l,mva_v]=irf_minvar(irf_tlim(E,tint));

    gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,mva_l,mva_v,f_highpass);      
    imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'eye/' irf_time(tint(1),'epoch2iso') '_' num2str(f_highpass) '_dir_.gif'],'DelayTime',0.01,'LoopCount',inf);
end
%% Plot correlation profiles
% cmaxs
%try fig = gcf; delete(fig.Children(:)); 1; end
clim = 0.8; 
ind_plot = find(cmaxs(:,2)>clim);
hca=subplot(2,1,1);
plot(hca,1:3:360,cs(ind_plot,:)')
hca.YLabel.String = ['C (C_{max> ' num2str(clim) '})'];
hca.XLabel.String = ['\theta'];
hca.XLim = [0 360];

hca=subplot(2,1,2); hold(hca,'on');

for kk = torow(ind_plot)
    % Find minimum
    indMin = find(cs(kk,:)==min(cs(kk,:)));
    plot(hca,(1:3:360)-indMin*3,cs(kk,:)')
    hca.YLabel.String = ['C (C_{max> ' num2str(clim) '})'];
    hca.XLabel.String = ['theta-theta@C_{min}'];
    disp(['index = ' num2str(kk) ' indMin = ' num2str(indMin)])
end
box(hca,'on');  
hold(hca,'off');  

    

%% 