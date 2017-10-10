cd /Users/Cecilia/Data/Cluster/20020330/

t1ov = toepoch([2002 03 30 12 50 00]);
t2ov = toepoch([2002 03 30 15 40 00]);
t1 = toepoch([2002 03 30 13 05 00]);
t2 = toepoch([2002 03 30 13 15 00]);
t1mp = toepoch([2002 03 30 13 11 00]);
t2mp = toepoch([2002 03 30 14 00 00]);
tint = [t1 t2];
%% Load data
% Data downloaded for shorter time interval
%caa_download(tint,'C?_CP_EFW_L2_E');
%caa_download(tint,'C?_CP_EFW_L2_E3D_GSE');
%caa_download(tint,'C?_CP_STA_CWF_GSE');
%caa_download(tint,'C?_CP_FGM_FULL');
doLoad=1;
if doLoad
    load matlab
else
    % Load data
    sclist=1:4;
    c_eval('diE?inert = c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
    c_eval('diE? = c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'',''mat'');',sclist);
    c_eval('diB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_ISR2'',''mat'');',sclist);
    c_eval('diB?fgm = c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);

    c_eval('gseE?inert = c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''mat'');',sclist);
    %c_eval('diE? = c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'',''mat'');',sclist);
    c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
    c_eval('gseB?fgm = c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
    %c_eval('gseB?fgm=irf_abs(gseB?fgm);',sclist)

    c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?sta);',sclist)
    c_eval('gseB?=irf_abs(gseB?);',sclist)
    c_eval('gsmB?=c_coord_trans(''GSE'',''GSM'',gseB?,''cl_id'',?);',sclist);

    c_eval('gseR? = c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);

    c_eval('[diE?,ang?]=irf_edb(diE?,diB?fgm);',sclist);
    c_eval('gseE?=c_coord_trans(''ISR2'',''GSE'',diE?,''cl_id'',?);',sclist);
    %gseE1=c_coord_trans('DSI','GSE',diE1,'cl_id',1);

    c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',sclist);
    c_eval('scpNe?=c_efw_scp2ne(P?);',sclist);
    c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',sclist);

    c_eval('parTe?=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
    %c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_MOMENTS'',''mat'');',3:4);

    units = irf_units;
    c_eval('flh? = gseB?(:,[1 5]);',sclist)
    Te=2000;
    c_eval('flh?(:,2) = flh?(:,2)*1e-9*units.e/sqrt(units.me*units.mp)/2/pi;',sclist)
    c_eval('vA? = irf_tappl(irf_plasma_calc(gseB?,scpNe?,0,Te,parTe?,''Va''),''*1e-3'');',sclist)
    c_eval('rhoe? = irf_plasma_calc(gseB?,scpNe?,0,Te,Te,''Roe'');',sclist)
    c_eval('Le? = irf_plasma_calc(gseB?,scpNe?,0,Te,Te,''Le'');',sclist)

    c_eval('beta? = irf_multiply(1,rhoe?,2,Le?,-1);',sclist)
end

%% Calculate wavelets and ftts
tzoom=toepoch([2002 03 30 13 11 30;2002 03 30 13 13 30])';
nfft = 512;
sfreq = 450;
overlap = 0.50;
fftE2 = irf_powerfft(gseE2,nfft,sfreq,[overlap]);
fftB2 = irf_powerfft(gseB2sta,nfft,sfreq,[overlap]);
wavE2 = irf_wavelet(irf_tlim(gseE2,tzoom),'Fs',sfreq);
c_eval('wavdiE? = irf_wavelet(irf_tlim(diE?(:,1:3),tzoom),''Fs'',sfreq);',2);
wavB2 = irf_wavelet(irf_tlim(gseB2sta,tzoom),'Fs',sfreq);
wavpolE2 = irf_wavelet(irf_tlim(irf_dot(gseE2(:,1:4),[-1 -0.06 0.04]),tzoom),'Fs',sfreq);
%% Plot data, Matlab_R2014a
h = irf_plot(8);
isub=1;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diB2fgm)
    ylabel(hca,'B_{ISR2,FGM} [mV/m]')  
end
if 1 % wave magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diB2sta,10,0))
    ylabel(hca,'B_{ISR2,STA} [mV/m]')  
end
if 1 % wave electric field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diE2,10,0))  
    ylabel(hca,'E_{ISR2} [mV/m]')  
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,fftE2,'log');
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,fftB2,'log');
    ylabel(hca,'f_B [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavE2,'log');
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavB2,'log');
    ylabel(hca,'f_B [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % ang between spin plane and B
    hca=h(isub); isub=isub+1;
    irf_plot(hca,[gseE2(:,1) ang2])
    ylabel(hca,'\theta_{sp-B}')
end
if 0 % wave electric field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diE2inert,10,0))
end

irf_zoom(h,'x',tzoom)
%% Plot data
h = irf_plot(7);
isub=1;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diB1fgm)
    ylabel(hca,'B_{ISR2,FGM} [mV/m]')  
end
if 1 % wave magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diB1sta,10,0))
    ylabel(hca,'B_{ISR2,STA} [mV/m]')  
end
if 1 % wave electric field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diE1,10,0))  
    ylabel(hca,'E_{ISR2} [mV/m]')  
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavdiE1.t,wavdiE1.p{1},wavdiE1.f);
    ylabel(hca,'f_{diEx} [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavdiE1.t,wavdiE1.p{2},wavdiE1.f);
    ylabel(hca,'f_{diEy} [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavpolE1.t,wavpolE1.p,wavpolE1.f);
    ylabel(hca,'f_{Epol} [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavE1.t,wavE1.p{1},wavE1.f);
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavE1.t,wavE1.p{2},wavE1.f);
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavE1.t,wavE1.p{3},wavE1.f);
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavB1,'log');
    ylabel(hca,'f_B [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end

irf_zoom(h,'x',tzoom)
%% Plot data, Matlab_R2014b
h = irf_plot(8);
isub=1;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diB1fgm)
    ylabel(hca,'B_{ISR2,FGM} [mV/m]')  
end
if 1 % wave magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diB1sta,10,0))
    ylabel(hca,'B_{ISR2,STA} [mV/m]')  
end
if 1 % wave electric field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diE1,10,0))  
    ylabel(hca,'E_{ISR2} [mV/m]')  
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,fftE1,'log');
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,fftB1,'log');
    ylabel(hca,'f_B [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavE1,'log');
    ylabel(hca,'f_E [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavB1,'log');
    ylabel(hca,'f_B [Hz]')
    hold(hca,'on')
    irf_plot(hca,wlh)
end
if 1 % ang between spin plane and B
    hca=h(isub); isub=isub+1;
    irf_plot(hca,[gseE1(:,1) ang1])
    ylabel(hca,'\theta_{sp-B}')
end
if 0 % wave electric field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_filt(diE1inert,10,0))
end

irf_zoom(h,'x',tzoom)

%% List of matching times, individual matching
t1 = toepoch([2002 03 30 13 11 44.00]);
t2 = toepoch([2002 03 30 13 11 44.30]);
%t1 = toepoch([2002 03 30 13 11 46.00]); t2 = toepoch([2002 03 30 13 11 46.10]);
t1 = toepoch([2002 03 30 13 11 51.50]); t2 = toepoch([2002 03 30 13 11 51.60]);
t1 = toepoch([2002 03 30 13 12 31.38]); t2 = toepoch([2002 03 30 13 12 31.46]);
%t1 = toepoch([2002 03 30 13 11 46.90]);
%t2 = toepoch([2002 03 30 13 11 47.10]);
t1 = toepoch([2002 03 30 13 11 50.82]); t2 = toepoch([2002 03 30 13 11 50.92]);
dt=0;
t1=t1-dt;
t2=t2-dt;
% Try matching
% Direction
angles=1:3:360;
f_highpass=50;
c_eval('E = irf_tlim(gseE?,[t1 t2]);',2)
c_eval('B = irf_tlim(gseB?,[t1 t2]);',2);
[x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
direction=x(i_dir,:);

% Velocity and density
n=linspace(0.01,0.1,100);
v=linspace(200,10000,30);
[corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
i_v=find(corr_v(:,1)==min(corr_v(:,1)));
velocity=v(i_v);

% Figures
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En);      
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,['mygif_dir' num2str(f_highpass) '_' irf_time(t1,'epoch2iso') '.gif'],'DelayTime',0.01,'LoopCount',inf);
%% Make gif of individual matching
i_n=50; % if more than one densitiy, choose one by specifying index
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);

figure; h=axes;
axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);


%% make minimum variance coordinate system
% mva from C2 i think
% from minimum variance analysis
v1 = [-0.1553 -0.2309 0.9605]; % max in GSE. 
v2 = -[-0.3610 0.9183 0.1623]; % inter in GSE. 
v3 = -[-0.9195 -0.3216 -0.2260]; % min in GSE. 

v1 = [-0.1654 -0.1976 0.9662]; % max in GSE. GSM: -0.1654   -0.4318    0.8866
v2 = -[-0.0774 0.9793 0.1870]; % inter in GSE. GSM: 0.0774   -0.9020   -0.4248
v3 = -[-0.9832 -0.0439 -0.1772]; % min in GSE. GSM: 0.0774   -0.9020   -0.4248
sclist=1:4;
c_eval('minR?=irf_lmn(gseR?,v1,v2,v3);',sclist);
c_eval('minE?=irf_lmn(gseE?,v1,v2,v3);',sclist);
c_eval('minB?=irf_lmn(gseB?(:,1:4),v1,v2,v3);',sclist);

v1 = [-0.1654 -0.1976 0.9662]; % max in GSE. GSM: -0.1654   -0.4318    0.8866
v2 = -[-0.0774 0.9793 0.1870]; % inter in GSE. GSM: 0.0774   -0.9020   -0.4248
v3 = -[-0.9832 -0.0439 -0.1772]; % min in GSE. GSM: 0.0774   -0.9020   -0.4248
sclist=1:4;
c_eval('slowB?=irf_add(1,gseB?,-1,irf_filt(gseB?,5,0));',sclist)
c_eval('facR?=irf_lmn(gseR?,slowB?,v3,''L'');',sclist);
c_eval('facE?=irf_lmn(gseE?,slowB?,v3,''L'');',sclist);
c_eval('facB?=irf_lmn(gseB?,slowB?,v3,''L'');',sclist);

%% run irf_ebsp, not complete, didn?t even start
tdisc = toepoch([2002 03 30 13 11 42;2002 03 30 13 11 49])'+[-20 + 20];
ebsp = irf_ebsp(irf_tlim(minE2,tdisc),irf_filt(irf_tlim(minB2,tdisc),10,0),irf_tlim(minB2,tdisc),irf_tlim(minB2,tdisc),[1 0 0],'pc12','fac')

%% estimate current sheet speed with [V, dV] = c_4_v_xcorr(tint,B1,B2,B3,B4,R1,R2,R3,R4)
tdisc = toepoch([2002 03 30 13 11 42;2002 03 30 13 11 49])';
step=6;
[V, dV] = c_4_v_xcorr(tdisc,gseB1(1:step:end,:),gseB2(1:step:end,:),gseB3(1:step:end,:),gseB4(1:step:end,:),...
                            gseR1(1:step:end,:),gseR2(1:step:end,:),gseR3(1:step:end,:),gseR4(1:step:end,:));
[V, dV] = c_4_v_xcorr(tdisc,minB1(1:step:end,:),minB2(1:step:end,:),minB3(1:step:end,:),minB4(1:step:end,:),...
                            minR1(1:step:end,:),minR2(1:step:end,:),minR3(1:step:end,:),minR4(1:step:end,:));

%c_eval('dt?! = (cn_average2(irf_tlim(minR!,tdisc),1)-cn_average2(irf_tlim(minR?,tdisc),1))./V;',1:4,1);
%dt = [dt11 dt21 dt31 dt41]

%from gui
dt = [0.00  1.65  2.96 -0.45]; % dt=[0.00  1.65  2.96 -0.45]
V = 31.2*[ -0.94 -0.17 -0.31]; % km/s GSE, strange, does not at all fit to the c_4_v_xcorr method

dt = [0.00      1.57      2.47     -0.64];
V = 29.3 * [-0.93 -0.35 -0.10];

minV=irf_lmn(V,v1,v2,v3); 

% make N position vector of discontinuity.
% choose one time when C1 exactly crosses the discontinuity
t0 = toepoch([2002 03 30 13 11 44.75]);
Rcs = minR1;
Rcs = [Rcs(:,1) Rcs(:,[2:4])+(Rcs(:,1)-t0)*minV];
%Rrel = (Rcs(:,1)-t0)*minV; % km
R0 = cn_average2(irf_tlim(minR1(:,:),t0+[-0.05 0.05]),3);
minRcs = minR1(:,[1 4]);
%% Calculate wavelets and ftts with LMN system
tzoom=toepoch([2002 03 30 13 11 30;2002 03 30 13 13 30])';
nfft = 512;
sfreq = 450;
overlap = 0.50;
sc = 1:2;
c_eval('wavLE? = irf_wavelet(irf_tlim(minE?(:,[1 2]),tzoom),''Fs'',sfreq);',sc);
c_eval('wavME? = irf_wavelet(irf_tlim(minE?(:,[1 3]),tzoom),''Fs'',sfreq);',sc);
c_eval('wavNE? = irf_wavelet(irf_tlim(minE?(:,[1 4]),tzoom),''Fs'',sfreq);',sc);
c_eval('wavLB? = irf_wavelet(irf_tlim(minB?(:,[1 2]),tzoom),''Fs'',sfreq);',sc);
c_eval('wavMB? = irf_wavelet(irf_tlim(minB?(:,[1 3]),tzoom),''Fs'',sfreq);',sc);
c_eval('wavNB? = irf_wavelet(irf_tlim(minB?(:,[1 4]),tzoom),''Fs'',sfreq);',sc);

%% plot and compare E and B wavelets
h = irf_plot(7);
isub=1;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,minB2)
    ylabel(hca,'B2 [nT]')  
    irf_legend(hca,{'L','M','N'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLE2,'log');
    ylabel(hca,'f_{E2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'L'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLB2,'log');
    ylabel(hca,'f_{B2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'L'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavME2,'log');
    ylabel(hca,'f_{E2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'M'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavMB2,'log');
    ylabel(hca,'f_{B2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'M'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavNE2,'log');
    ylabel(hca,'f_{E2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'N'},[0.95 0.95])
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavNB2,'log');
    ylabel(hca,'f_{B2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'N'},[0.95 0.95])
end

if 0 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,minB3)
    ylabel(hca,'B3 [nT]')  
    irf_legend(hca,{'L','M','N'},[0.95 0.95])
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLE3,'log');
    ylabel(hca,'f_{E3} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh3)
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLB3,'log');
    ylabel(hca,'f_{B3} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh3)
end

if 0 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,minB4)
    ylabel(hca,'B4 [nT]')  
    irf_legend(hca,{'L','M','N'},[0.95 0.95])
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLE4,'log');
    ylabel(hca,'f_{E4} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh4)
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLB4,'log');
    ylabel(hca,'f_{B4} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh4)
end
for kk = [2 4 6]; % E   
    %ch = colorbar('peer',h(kk));
    set(h(kk),'clim',[-5 1.4])    
end
for kk = [3 5 7]; % B
    %ch = colorbar('peer',h(kk));
    set(h(kk),'clim',[-5 1.8])    
end
irf_zoom(h,'x',tzoom)
irf_plot_axis_align
%% moving matching correlation!
ta1 = toepoch([2002 03 30 13 11 40.00]); % start time of wave packet
ta2 = toepoch([2002 03 30 13 11 52.00]); % end time of wave packet
%ta1 = toepoch([2002 03 30 13 12 25.00]); % start time of wave packet
%ta2 = toepoch([2002 03 30 13 12 40.00]); % end time of wave packet
%ta1 = toepoch([2002 03 30 13 13 00.00]); % start time of wave packet
%ta2 = toepoch([2002 03 30 13 13 22.00]); % end time of wave packet
%ta1 = toepoch([2002 03 30 13 11 40.00]); % start time of wave packet
tint_a = [ta1 ta2];

% filtering done att 0.8*f_lh
f_lim = 0.8;
sclist=1:4;
c_eval('f_filt? = flh?; f_filt?(:,2) = f_filt?(:,2)*f_lim;',sclist);

% number of wave periods to include in each moving bin, try 5?
% just take it for an average flh.. flh = 60 gives T = nTlh*1/flh
nTlh = 5;
flh_av = 35;
T = nTlh/flh_av;
%nT = diff(tint_a)/T;
t0s = ta1 + 0:0.5*T:ta2;
nT = numel(t0s);
angles=1:3:360;

c_eval('minBac? = irf_filt(minB?(:,1:4),flh_av,0);',sclist);

sclist = 2;

for sc = sclist;
    % initialize variables
    c_eval('ts? = zeros(nT,1);',sc);
    c_eval('ks? = zeros(nT,4);',sc);    
    c_eval('ns? = zeros(nT,4);',sc);
    c_eval('Bs? = zeros(nT,4);',sc);
    c_eval('cs? = zeros(nT,numel(angles));',sc);
    c_eval('cmaxs? = zeros(nT,2);',sc);
    c_eval('mB? = zeros(nT,2);',sc);
    c_eval('mE? = zeros(nT,2);',sc);
    c_eval('mintEdt? = zeros(nT,2);',sc);
    c_eval('B0s? = zeros(nT,2);',sc);
    c_eval('Ts? = zeros(nT,2);',sc);
    c_eval('tints? = zeros(nT,3);',sc);
    
    c_eval('kmvas? = zeros(nT,4);',sc);
    c_eval('mvaang? = zeros(nT,2);',sc);
    c_eval('l12_? = zeros(nT,2);',sc);
    c_eval('l23_? = zeros(nT,2);',sc);    
    %f_filt = zeros(nT-1,2);

    %sc = 2;
    for kk = 1:nT  
        %tint_loc = [tints(kk) tints(kk)+T];

        c_eval('flh_loc = cn_toepoch(t0s(kk),flh?);',sc);
        T_loc = nTlh/flh_loc(2);
        %T_loc = nTlh/flh_av;
        tint_loc = t0s(kk) + 0.5*T_loc*[-1 1];    
    
        %f_filt(kk,:) = cn_average2(irf_tlim(f_lh,tint_loc),3);
        %t_mean = sum(tint_loc)*0.5;
        t_mean = t0s(kk);
        c_eval('Ts?(kk,:) = [t_mean T_loc];',sc);
        c_eval('tints?(kk,:) = [t_mean tint_loc];',sc);

        % take out lower hybrid frequency for filtering.
        c_eval('f_filt_loc = cn_average2(irf_tlim(f_filt?,tint_loc),1);',sc);
        c_eval('E = irf_tlim(minE?,tint_loc);',sc)
        c_eval('B = irf_tlim(minB?(:,1:4),tint_loc);',sc);
                
        [x,y,z,corr_dir,intEdt,Bz,B0,dEk,dEn,Ek,En]=irf_match_phibe_dir(B,E,angles,f_filt_loc);
        c_eval('mB?(kk,:) = [t_mean mean(abs(Bz(:,2)))];',sc);
        
        c_eval('cmaxs?(kk,:) = [t_mean max(corr_dir(:,1))];',sc);
        i_dir = find(corr_dir(:,1)==max(corr_dir(:,1)));    
        c_eval('mE?(kk,:) = [t_mean mean(abs(dEk(2:end,i_dir+1)))];',sc);
        c_eval('mintEdt?(kk,:) = [t_mean mean(abs(intEdt(2:end,i_dir+1)))];',sc);

        c_eval('ks?(kk,:) = [t_mean x(i_dir,:)];',sc);
        c_eval('ns?(kk,:) = [t_mean y(i_dir,:)];',sc);
        c_eval('Bs?(kk,:) = [t_mean z(i_dir,:)];',sc);

        c_eval('B0s?(kk,:) = [t_mean B0];',sc);
        c_eval('cs?(kk,:) = corr_dir'';',sc);
        
        c_eval('[~,l,v]=irf_minvar(irf_tlim(minBac?,tint_loc));',sc);
        c_eval('mvaang?(kk,:) = [t_mean (abs((z(i_dir,:)*v(3,:)'')))];',sc);
        c_eval('kmvas?(kk,:) = [t_mean torow(v(3,:))];',sc);
        c_eval('l12_?(kk,:) = [t_mean l(1)/l(2)];',sc);
        c_eval('l23_?(kk,:) = [t_mean l(2)/l(3)];',sc);
        %ang2B_vec(yind)=acosd(abs((z(i_dir,:)*v(3,:)')));
        %ang2x_vec(yind)=acosd(abs(([1 0 0]*v(3,:)')));
    end

    c_eval('EB? = irf_multiply(1e-3/1e-9,mE?,1,mB?,-1);',sc);

    c_eval('scpmNe?=irf_resamp(scpNe?,ks?(:,1));',sc);
    c_eval('mv? = [ks?(:,1) 1e-3*B0s?(:,2)*1e-9.*mB?(:,2)*1e-9./(units.e*units.mu0.*scpmNe?(:,2)*1e6.*mintEdt?(:,2)*1e-3)];',sc);
end
%% plot results, many panels 
h=irf_plot(10);
isub=1;
coord_leg = {'L','M','N'};
sc = 2;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_abs(minB?))',sc);
    ylabel(hca,'B [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 0 % electron density
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{peaNe?,scpNe?},''comp'')',sc);
    ylabel(hca,'n_e [cc]')
    irf_legend(hca,{'PEA','EFW'},[1 0.95])
end
if 0 % DPFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dEF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    %caxis([5.9 8.6]);
    set(hca,'yscale','log','ylim',[11 3e4]);
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
end
if 1 % electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,minE?)',sc);
    ylabel(hca,'E [mV/m]') 
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 1 % wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_filt(minB?,30,0))',sc);
    ylabel(hca,'\delta B [nT]') 
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 0 % electron temperature
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,parTe?)',sc);
    ylabel(hca,'T_e [eV]')
    %irf_legend(hca,{'PEA','EFW'},[1 0.95])
end
if 0 % electron plasma beta
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,beta?)',sc);
    set(hca,'yscale','log','ytick',[ 1e0 1e1 1e2 1e3 1e4],'ylim',[1e2 5e4])    
    ylabel(hca,'\beta ')
    %irf_legend(hca,{'PEA','EFW'},[1 0.95])
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavLB2,'log');
    ylabel(hca,'f_{B2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'L'},[0.95 0.95])
end
if 0 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,wavME2,'log');
    ylabel(hca,'f_{E2} [Hz]')
    hold(hca,'on')
    irf_plot(hca,flh2)
    irf_legend(hca,{'M'},[0.95 0.95])
end

if 1 % k direction
    hca=h(isub); isub=isub+1;
    corr_lim=0.0;
    c_eval('plot_k = ks?;',sc);
    c_eval('plot_k(cmaxs?(:,2)<corr_lim,:)=NaN;',sc);
    irf_plot(hca,plot_k)
    ylabel(hca,'k') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 1 % normal direction
    hca=h(isub); isub=isub+1;
    c_eval('plot_n = ns?;',sc); 
    c_eval('plot_n(cmaxs?(:,2)<corr_lim,:)=NaN;',sc);
    irf_plot(hca,plot_n)
    ylabel(hca,'n') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 0 % B0 direction
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,Bs?)',sc);
    ylabel(hca,'B') 
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,ks?(:,1),cs?'',angles);',sc);
    set(hca,'clim',[-1 1])
    ch = colorbar('peer',hca);
    ylabel(ch,'correlation')
    ylabel(hca,{'\theta to','(B\times[0 0 1])\times B'})
end
if 1 % max correlation
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,cmaxs?);',sc);
    ylabel(hca,{'corr_{max}'})
    irf_zoom(hca,'y',[0 1])
end
if 0 % <|E|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mE?);',sc);
    ylabel(hca,{'<|E|>','[mV/m]'})
end
if 0 % <|Phi|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mintEdt?);',sc);
    ylabel(hca,{'<|Edt|>','[(mV/m)s]'})
end
if 0 % <|B|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mB?);',sc);
    ylabel(hca,{'<|B|>','[nT]'})
end
if 0 % vA
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,vA?);',sc);
    ylabel(hca,{'v_A','[km/s]'})
    %irf_zoom(hca,'y',[0 50000])
end
if 0 % E/B ratio
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,irf_tappl(EB?,''*1e-3''));',sc);
    ylabel(hca,{'<|E|>/<|B|>','[km/s]'})
    %irf_zoom(hca,'y',[0 50000])
end

if 1 % <|v|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mv?);',sc);
    ylabel(hca,{'<|v_{EH}|>','[km/s]'})
    irf_zoom(hca,'y',[0 1000])
end
if 1 % <|v|> /flh
    hca=h(isub); isub=isub+1;
    %c_eval('irf_plot(hca,mv?);',sc);
    %c_eval('irf_plot(hca,irf_multiply(1e3,irf_multiply(1,mv?,1,flh?,-1),1,rhoe?,-1));',sc);
    c_eval('irf_plot(hca,{irf_multiply(1,mv?,1,flh?,-1),irf_tappl(rhoe?,''*1e-3'')},''comp'');',sc);
    ylabel(hca,{'<|v_{EH}|>/f_{LH}','[km]'})
    irf_zoom(hca,'y',[0 30])
end
if 0 % E/B/v_A ratio
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,irf_multiply(1e-3,EB?,1,vA?,-1));',sc);
    ylabel(hca,{'<|E|>/<|B|>/','v_A'})
    %irf_zoom(hca,'y',[0 50000])
end
if 1 % E/B/v_EH ratio
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,irf_multiply(1e-3,EB?,1,irf_tappl(mv?,''*1e0''),-1));',sc);
    ylabel(hca,{'<|E|>/<|B|>/','v_{EH}'})
    irf_zoom(hca,'y',[0 400])
end
if 0 % ang between spin plane and B
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,[gseE?(:,1) ang?])',sc);
    ylabel(hca,'\theta_{sp-B}')
end
if 0 % flh and ffilt
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{flh?,f_filt?},''comp'')',sc);
    ylabel(hca,'f [Hz]')
    irf_legend(hca,{'f_{LH}','f_{filt}'},[0.95 0.9])
end
if 0 % local time period
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,Ts?)',sc);
    ylabel(hca,'T [s]')
    %irf_legend(hca,{'f_{LH}','f_{filt}'},[0.95 0.9])
end

irf_zoom(h,'x',tint_a+[+3 0])
irf_plot_axis_align
%% correlation comparison between the three sc
h=irf_plot(6);
isub = 1;
coord_leg = {'L','M','N'};
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(minB2))
    ylabel(hca,'B_2 [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 1
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,ks2(:,1),cs2',angles);    
    set(hca,'clim',[-1 1]);
    ch = colorbar('peer',hca);
    ylabel(hca,'\theta_2')
end

if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(minB3))
    ylabel(hca,'B_3 [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 1
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,ks3(:,1),cs3',angles);
    set(hca,'clim',[-1 1]);
    ch = colorbar('peer',hca);
    ylabel(hca,'\theta_3')
end

if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(minB4))
    ylabel(hca,'B_4 [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 1
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,ks4(:,1),cs4',angles);
    set(hca,'clim',[-1 1]);
    ch = colorbar('peer',hca);
    ylabel(hca,'\theta_4')
end

irf_zoom(h,'x',tint_a+[3 0])
irf_plot_axis_align
%% debug v_eh, fine now
hhh=irf_plot({B0s,mB,ns,Phi_mean,mv,irf_multiply(1,Phi_mean,1,mv,1)});
ylabel(hhh(1),'B [nT]')
ylabel(hhh(2),'\delta B [nT]')
ylabel(hhh(3),'n [cc]')
ylabel(hhh(4),'\int Edt [mV m^{-1} s]')
ylabel(hhh(5),'v [km/s]')
ylabel(hhh(6),'\int Edt \times v [V]')
%% make gif to show how k,n,B changes in LMN-system, 1 sc
% x = L
% y = M
% z = N
% make patch to show LM-plane
xpl = 2*[-1 -1 1 1];
ypl = 2*[1 -1 -1 1];
zpl = 0*[-1 -1 1 1];

% to print with opacity, figure renderer must be set to OpenGl
nPanels = 4;
for kk = 1:nPanels; h(kk) = subplot(nPanels,1,kk); end
isub = 1;

% set up plot for quivers
hca = h(isub); isub=isub+1;
axes(h(1)); % set current axes
patch(xpl,ypl,zpl);
hp = findobj(gcf,'type','patch');
set(hp,'facealpha',0.05)
xlabel(hca,'L')
ylabel(hca,'M')
zlabel(hca,'N')
hold(hca,'on')
view(hca,[1 1 0.3])
nt = numel(ks(:,1));
xlabelspatch = get(hca,'xticklabel');
xtickspatch = get(hca,'xtick');
box(hca,'on')

% set up plot for time series data
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(minB2))
    ylabel(hca,'B [nT]')
    irf_legend(hca,{'L','M','N'},[1 0.95])
end

if 1 % field spectrogram
    hca=h(isub); isub=isub+1;
    irf_spectrogram(hca,ks(:,1),cs',angles);
    ylabel(hca,'angles')
end
if 1 % max correlation
    hca=h(isub); isub=isub+1;
    irf_plot(hca,cmaxs);
    ylabel(hca,'corr_{max}')
    irf_zoom(hca,'y',[0 1])
end

irf_zoom(h(2:4),'x',tint_a+[3 0])
set(h(1),'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'ztick',[],'zticklabel',[])
%% continuation
for kk = 1:nt    
    lines=findobj(gcf,'type','line');
    delete(lines);
    quiver3(hca,0,0,0,ks(kk,2),ks(kk,3),ks(kk,4),'k')
    quiver3(hca,0,0,0,ns(kk,2),ns(kk,3),ns(kk,4),'b')
    quiver3(hca,0,0,0,Bs(kk,2),Bs(kk,3),Bs(kk,4),'r')
    pause(0.05)
    irf_plot([ks(kk,1)])
end

%% plot three sc at the same time
% have min R?
% make locR?
minR0 = [minR2(:,1) (minR2(:,2:4)+minR3(:,2:4)+minR4(:,2:4))/3];
%R0 = irf_tlim(R0,tint_a);
c_eval('locR? = irf_add(1,minR?,-1,minR0);',2:4);
c_eval('R? = cn.mean(irf_tlim(locR?,tint_a),1);',2:4);

% x = L
% y = M
% z = N
% make patch to show LM-plane
xmax = 1.5*max(abs([R2(1) R3(1) R4(1)])); 
ymax = 1.5*max(abs([R2(2) R3(2) R4(2)]));
zmax = 1.5*max(abs([R2(3) R3(3) R4(3)]));
xpl = xmax*[-1 -1 1 1];
ypl = ymax*[1 -1 -1 1];
zpl = 0*[-1 -1 1 1];

nPanels = 1;
for kk = 1:nPanels; h(kk) = subplot(nPanels,1,kk); end
isub = 1;
hca = h(isub); isub=isub+1;

% plot spacecraft positions
plot3(hca,R2(1),R2(2),R2(3),'rs',R3(1),R3(2),R3(3),'gs',R4(1),R4(2),R4(3),'bs','markersize',20)

% make patch to vizualize plane of current sheet, but not necessarily
% absolute placement in the normal direction
% to print with opacity, figure renderer must be set to OpenGl
axes(h(1)); % set current axes
patch(xpl,ypl,zpl);
hp = findobj(gcf,'type','patch');
set(hp,'facealpha',0.05)
xlabel(hca,'L')
ylabel(hca,'M')
zlabel(hca,'N')
hold(hca,'on')
view(hca,[1 1 0.3])
%view(hca,[0 1 0])
axis equal
set(hca,'xlim',xmax*[-1 1],'ylim',ymax*[-1 1],'zlim',zmax*[-1 1]);
xlabelspatch = get(hca,'xticklabel');
xtickspatch = get(hca,'xtick');
box(hca,'on')
%
hold(gca,'on')
nt = numel(ks2(:,1));
qsize = 20;
for kk = 1:nt    
    title(irf_time(ks2(kk),'epoch2isoshort'))
    lines=findobj(gcf,'type','line');
    delete(lines);
    quiver3(hca,R2(1),R2(2),R2(3),ks2(kk,2),ks2(kk,3),ks2(kk,4),qsize,'k')
    quiver3(hca,R2(1),R2(2),R2(3),ns2(kk,2),ns2(kk,3),ns2(kk,4),qsize,'b')
    quiver3(hca,R2(1),R2(2),R2(3),Bs2(kk,2),Bs2(kk,3),Bs2(kk,4),qsize,'r')
    quiver3(hca,R3(1),R3(2),R3(3),ks3(kk,2),ks3(kk,3),ks3(kk,4),qsize,'k')
    quiver3(hca,R3(1),R3(2),R3(3),ns3(kk,2),ns3(kk,3),ns3(kk,4),qsize,'b')
    quiver3(hca,R3(1),R3(2),R3(3),Bs3(kk,2),Bs3(kk,3),Bs3(kk,4),qsize,'r')
    quiver3(hca,R4(1),R4(2),R4(3),ks4(kk,2),ks4(kk,3),ks4(kk,4),qsize,'k')
    quiver3(hca,R4(1),R4(2),R4(3),ns4(kk,2),ns4(kk,3),ns4(kk,4),qsize,'b')
    quiver3(hca,R4(1),R4(2),R4(3),Bs4(kk,2),Bs4(kk,3),Bs4(kk,4),qsize,'r')
    plot3(hca,R2(1),R2(2),R2(3),'rs',R3(1),R3(2),R3(3),'gs',R4(1),R4(2),R4(3),'bs','markersize',20)
    pause(0.1)
    %irf_plot([ks(kk,1)])
end

%%
% make N position vector of discontinuity.
% choose one time when C1 exactly crosses the discontinuity
t0 = toepoch([2002 03 30 13 11 45.2]);
Rcs = minR1;
Rcs = [Rcs(:,1) Rcs(:,[2:4])+(Rcs(:,1)-t0)*minV];
%Rrel = (Rcs(:,1)-t0)*minV; % km
R0 = cn_average2(irf_tlim(minR1(:,:),t0+[-0.05 0.05]),3);
minRcs = minR1(:,[1 4]);

%% plot three sc at the same time + moving current sheet
% have min R?
% make locR?
minR0 = [minR2(:,1) (minR1(:,2:4)+minR2(:,2:4)+minR3(:,2:4)+minR4(:,2:4))/4];
%R0 = irf_tlim(R0,tint_a);
c_eval('locR? = irf_add(1,minR?,-1,minR0);',1:4);
locRcs = irf_add(1,Rcs,-1,minR0); locRcs = irf_resamp(locRcs,ks2(:,1));
c_eval('R? = cn.mean(irf_tlim(locR?,tint_a),1);',1:4);

% x = L
% y = M
% z = N
% make patch to show LM-plane
xmax = 2*max(abs([R1(1) R2(1) R3(1) R4(1)])); 
ymax = 2*max(abs([R1(2) R2(2) R3(2) R4(2)]));
zmax = 1*max(abs([R1(3) R2(3) R3(3) R4(3) locRcs(:,4)']));
xpl = xmax*[-1 -1 1 1];
ypl = ymax*[1 -1 -1 1];


nPanels = 1;
for kk = 1:nPanels; h(kk) = subplot(nPanels,1,kk); end
isub = 1;
hca = h(isub); isub=isub+1;

% plot spacecraft positions
plot3(hca,R1(1),R1(2),R1(3),'ks',R2(1),R2(2),R2(3),'rs',R3(1),R3(2),R3(3),'gs',R4(1),R4(2),R4(3),'bs','markersize',20)

%
%
box(hca,'on')
hold(gca,'on')
nt = numel(ks2(:,1));
qsize = 20;
for kk = 40:10:160;nt  ;  
    hp = findobj(gcf,'type','patch'); delete(hp);
    % make patch to vizualize plane of current sheet, but not necessarily
    % absolute placement in the normal direction
    % to print with opacity, figure renderer must be set to OpenGl
    zpl = locRcs(kk,4)*[1 1 1 1];
    axes(h(1)); % set current axes
    patch(xpl,ypl,zpl,[0 1 0]);
    hp = findobj(gcf,'type','patch');
    set(hp,'facealpha',0.05)
    xlabel(hca,'L')
    ylabel(hca,'M')
    zlabel(hca,'N')
    hold(hca,'on')
    view(hca,[1 1 0.7])
    %view(hca,[0 1 0])
    axis equal
    set(hca,'xlim',xmax*[-1 1],'ylim',ymax*[-1 1],'zlim',zmax*[-1 1]);
    xlabelspatch = get(hca,'xticklabel');
    xtickspatch = get(hca,'xtick');

    title({irf_time(ks2(kk),'epoch2isoshort'),['c_1=' num2str(cmaxs1(kk,2)) ',  c_2=' num2str(cmaxs2(kk,2)) ',  c_3=' num2str(cmaxs3(kk,2)) ',  c_4=' num2str(cmaxs4(kk,2))]})
    lines=findobj(gcf,'type','line');
    delete(lines);
    %quiver3(hca,R1(1),R1(2),R1(3),ks1(kk,2),ks1(kk,3),ks1(kk,4),qsize,'k')
    %quiver3(hca,R1(1),R1(2),R1(3),ns1(kk,2),ns1(kk,3),ns1(kk,4),qsize,'b')
    quiver3(hca,R1(1),R1(2),R1(3),Bs1(kk,2),Bs1(kk,3),Bs1(kk,4),qsize,'r','linewidth',2)
    quiver3(hca,R2(1),R2(2),R2(3),ks2(kk,2),ks2(kk,3),ks2(kk,4),qsize,'g','linewidth',2)
    quiver3(hca,R2(1),R2(2),R2(3),ns2(kk,2),ns2(kk,3),ns2(kk,4),qsize,'b','linewidth',2)
    quiver3(hca,R2(1),R2(2),R2(3),Bs2(kk,2),Bs2(kk,3),Bs2(kk,4),qsize,'r','linewidth',2)
    quiver3(hca,R3(1),R3(2),R3(3),ks3(kk,2),ks3(kk,3),ks3(kk,4),qsize,'g','linewidth',2)
    quiver3(hca,R3(1),R3(2),R3(3),ns3(kk,2),ns3(kk,3),ns3(kk,4),qsize,'b','linewidth',2)
    quiver3(hca,R3(1),R3(2),R3(3),Bs3(kk,2),Bs3(kk,3),Bs3(kk,4),qsize,'r','linewidth',2)
    quiver3(hca,R4(1),R4(2),R4(3),ks4(kk,2),ks4(kk,3),ks4(kk,4),qsize,'g','linewidth',2)
    quiver3(hca,R4(1),R4(2),R4(3),ns4(kk,2),ns4(kk,3),ns4(kk,4),qsize,'b','linewidth',2)
    quiver3(hca,R4(1),R4(2),R4(3),Bs4(kk,2),Bs4(kk,3),Bs4(kk,4),qsize,'r','linewidth',2)
    
    % projection on LM plane    
    quiver3(hca,R1(1),R1(2),-zmax,Bs1(kk,2),Bs1(kk,3),Bs1(kk,4)*0,qsize,'r')
    quiver3(hca,R2(1),R2(2),-zmax,ks2(kk,2),ks2(kk,3),ks2(kk,4)*0,qsize,'g')
    quiver3(hca,R2(1),R2(2),-zmax,ns2(kk,2),ns2(kk,3),ns2(kk,4)*0,qsize,'b')
    quiver3(hca,R2(1),R2(2),-zmax,Bs2(kk,2),Bs2(kk,3),Bs2(kk,4)*0,qsize,'r')
    quiver3(hca,R3(1),R3(2),-zmax,ks3(kk,2),ks3(kk,3),ks3(kk,4)*0,qsize,'g')
    quiver3(hca,R3(1),R3(2),-zmax,ns3(kk,2),ns3(kk,3),ns3(kk,4)*0,qsize,'b')
    quiver3(hca,R3(1),R3(2),-zmax,Bs3(kk,2),Bs3(kk,3),Bs3(kk,4)*0,qsize,'r')
    quiver3(hca,R4(1),R4(2),-zmax,ks4(kk,2),ks4(kk,3),ks4(kk,4)*0,qsize,'g')
    quiver3(hca,R4(1),R4(2),-zmax,ns4(kk,2),ns4(kk,3),ns4(kk,4)*0,qsize,'b')
    quiver3(hca,R4(1),R4(2),-zmax,Bs4(kk,2),Bs4(kk,3),Bs4(kk,4)*0,qsize,'r')
    
    % projection on NL plane
    quiver3(hca,R1(1),-ymax,R1(3),Bs1(kk,2),Bs1(kk,3)*0,Bs1(kk,4),qsize,'r')
    quiver3(hca,R2(1),-ymax,R2(3),ks2(kk,2),ks2(kk,3)*0,ks2(kk,4),qsize,'g')
    quiver3(hca,R2(1),-ymax,R2(3),ns2(kk,2),ns2(kk,3)*0,ns2(kk,4),qsize,'b')
    quiver3(hca,R2(1),-ymax,R2(3),Bs2(kk,2),Bs2(kk,3)*0,Bs2(kk,4),qsize,'r')
    quiver3(hca,R3(1),-ymax,R3(3),ks3(kk,2),ks3(kk,3)*0,ks3(kk,4),qsize,'g')
    quiver3(hca,R3(1),-ymax,R3(3),ns3(kk,2),ns3(kk,3)*0,ns3(kk,4),qsize,'b')
    quiver3(hca,R3(1),-ymax,R3(3),Bs3(kk,2),Bs3(kk,3)*0,Bs3(kk,4),qsize,'r')
    quiver3(hca,R4(1),-ymax,R4(3),ks4(kk,2),ks4(kk,3)*0,ks4(kk,4),qsize,'g')
    quiver3(hca,R4(1),-ymax,R4(3),ns4(kk,2),ns4(kk,3)*0,ns4(kk,4),qsize,'b')
    quiver3(hca,R4(1),-ymax,R4(3),Bs4(kk,2),Bs4(kk,3)*0,Bs4(kk,4),qsize,'r')
    
    % projection on NM plane
    quiver3(hca,-xmax,R1(2),R1(3),Bs1(kk,2)*0,Bs1(kk,3),Bs1(kk,4),qsize,'r')
    quiver3(hca,-xmax,R2(2),R2(3),ks2(kk,2)*0,ks2(kk,3),ks2(kk,4),qsize,'g')
    quiver3(hca,-xmax,R2(2),R2(3),ns2(kk,2)*0,ns2(kk,3),ns2(kk,4),qsize,'b')
    quiver3(hca,-xmax,R2(2),R2(3),Bs2(kk,2)*0,Bs2(kk,3),Bs2(kk,4),qsize,'r')
    quiver3(hca,-xmax,R3(2),R3(3),ks3(kk,2)*0,ks3(kk,3),ks3(kk,4),qsize,'g')
    quiver3(hca,-xmax,R3(2),R3(3),ns3(kk,2)*0,ns3(kk,3),ns3(kk,4),qsize,'b')
    quiver3(hca,-xmax,R3(2),R3(3),Bs3(kk,2)*0,Bs3(kk,3),Bs3(kk,4),qsize,'r')
    quiver3(hca,-xmax,R4(2),R4(3),ks4(kk,2)*0,ks4(kk,3),ks4(kk,4),qsize,'g')
    quiver3(hca,-xmax,R4(2),R4(3),ns4(kk,2)*0,ns4(kk,3),ns4(kk,4),qsize,'b')
    quiver3(hca,-xmax,R4(2),R4(3),Bs4(kk,2)*0,Bs4(kk,3),Bs4(kk,4),qsize,'r')
    
    plot3(hca,R1(1),R1(2),R1(3),'ks',R2(1),R2(2),R2(3),'rs',R3(1),R3(2),R3(3),'gs',R4(1),R4(2),R4(3),'bs','markersize',20)
    pause(0.1)
    %irf_plot([ks(kk,1)])
end

%% check wave polarization
% check wave polarization
v1 = [-0.14 -0.04 0.99]; % max
v2 = [-0.19 -0.98 -0.07]; % inter
v3 = [0.97 -0.20 0.13]; % min

% cn_average2(irf_tlim(minB3,toepoch([2002 03 30 13 12 30;2002 03 30 13 12 31])'),1)/norm(cn_average2(irf_tlim(minB3,toepoch([2002 03 30 13 12 30;2002 03 30 13 12 31])'),1))
B2hat = [0.5466    0.8360   -0.0486];
B3hat = [0.5732    0.8191    0.0249];
acosd(v3*B2hat') % angle between B and minBwave ("k" from k dot B=0)

% compare amplitudes of parBwave and perpBwave
ffilt=5;
B3wave = irf_filt(irf_tlim(minB3,toepoch([2002 03 30 13 12 30;2002 03 30 13 12 31])'),ffilt,0);

lmnB3wave = irf_lmn(B3wave,B3hat,[0 0 1],'L');
h=irf_plot(lmnB3wave);
irf_legend({'L','M','N'},[1 0.95])
%% compare tool with bool
% tool gives ks?
% bool gives kmvas?

h=irf_plot(11);
isub=1;
coord_leg = {'L','M','N'};
sc = 2;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_abs(minB?))',sc);
    ylabel(hca,'B [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 1 % electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,minE?)',sc);
    ylabel(hca,'E [mV/m]') 
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 1 % wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_filt(minB?,30,0))',sc);
    ylabel(hca,'\delta B [nT]') 
    irf_legend(hca,coord_leg,[1 0.95]) 
end
if 1 % k directions from bool and tool, 3 panels
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 2]),kmvas?(:,[1 2])},''comp'')',sc);
    ylabel(hca,'k_L') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 3]),kmvas?(:,[1 3])},''comp'')',sc);
    ylabel(hca,'k_M') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 4]),kmvas?(:,[1 4])},''comp'')',sc);
    ylabel(hca,'k_N') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
end
if 1 % angle between ktool and kbool
    hca=h(isub); isub=isub+1;
    c_eval('kkdot = irf_dot(ks?,kmvas?);',sc)
    anngg = [kkdot(:,1) acosd(kkdot(:,2))];
    irf_plot(hca,anngg);
    ylabel(hca,{'\theta_{kk}'})
    irf_zoom(hca,'y',[0 180])
    set(hca,'ytick',[0 30 60 90 120 150 180])
end
if 1 % angle between ktool and kbool
    hca=h(isub); isub=isub+1;
    %c_eval('kkdot = irf_dot(ks?,kmvas?);',sc)
    %anngg = [kkdot(:,1) acosd(kkdot(:,2))];
    c_eval('irf_plot(hca,mvaang?);',sc)
    ylabel(hca,{'\theta_{bool,B}'})
    irf_zoom(hca,'y',[0 180])
    set(hca,'ytick',[0 30 60 90 120 150 180])
end
if 1 % max correlation
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,cmaxs?);',sc);
    ylabel(hca,{'corr_{max}'})
    irf_zoom(hca,'y',[0 1])
end
if 1 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l12_?);',sc);
    ylabel(hca,{'L_{1}/L_{2}'})
    %irf_zoom(hca,'y',[0 1])
end
if 1 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l23_?);',sc);
    ylabel(hca,{'L_{2}/L_{3}'})
    %irf_zoom(hca,'y',[0 1])
end
if 0 % <|B|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mB?);',sc);
    ylabel(hca,{'<|B|>','[nT]'})
end
if 0 % E/B ratio
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,irf_tappl(EB?,''*1e-3''));',sc);
    ylabel(hca,{'<|E|>/<|B|>','[km/s]'})
    %irf_zoom(hca,'y',[0 50000])
end

irf_zoom(h,'x',tint_a+[+3 0])
irf_plot_axis_align
%% compare the EB polarizations
% tool gives ks?
% bool gives kmvas?

h=irf_plot(5);
isub=1;
coord_leg = {'L','M','N'};
sc = 2;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_abs(minB?))',sc);
    ylabel(hca,'B [nT]')
    irf_legend(hca,coord_leg,[1 0.95])
end
if 1 % MVA electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,minE?)',sc);
    ylabel(hca,'E [mV/m]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end
if 1 % MVA wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_filt(minB?,30,0))',sc);
    ylabel(hca,'\delta B [nT]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end
if 1 % FAC electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,facE?)',sc);
    ylabel(hca,'E_{fac} [mV/m]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end
if 1 % FAC wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_filt(facB?(:,1:4),30,0))',sc);
    ylabel(hca,'\delta B_{fac} [nT]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end

if 0 % k directions from bool and tool, 3 panels
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 2]),kmvas?(:,[1 2])},''comp'')',sc);
    ylabel(hca,'k_L') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 3]),kmvas?(:,[1 3])},''comp'')',sc);
    ylabel(hca,'k_M') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 4]),kmvas?(:,[1 4])},''comp'')',sc);
    ylabel(hca,'k_N') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
end
if 0 % angle between ktool and kbool
    hca=h(isub); isub=isub+1;
    c_eval('kkdot = irf_dot(ks?,kmvas?);',sc)
    anngg = [kkdot(:,1) acosd(kkdot(:,2))];
    irf_plot(hca,anngg);
    ylabel(hca,{'\theta_{kk}'})
    irf_zoom(hca,'y',[0 180])
    set(hca,'ytick',[0 30 60 90 120 150 180])
end
if 0 % angle between ktool and kbool
    hca=h(isub); isub=isub+1;
    %c_eval('kkdot = irf_dot(ks?,kmvas?);',sc)
    %anngg = [kkdot(:,1) acosd(kkdot(:,2))];
    c_eval('irf_plot(hca,mvaang?);',sc)
    ylabel(hca,{'\theta_{bool,B}'})
    irf_zoom(hca,'y',[0 180])
    set(hca,'ytick',[0 30 60 90 120 150 180])
end
if 0 % max correlation
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,cmaxs?);',sc);
    ylabel(hca,{'corr_{max}'})
    irf_zoom(hca,'y',[0 1])
end
if 0 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l12_?);',sc);
    ylabel(hca,{'L_{1}/L_{2}'})
    %irf_zoom(hca,'y',[0 1])
end
if 0 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l23_?);',sc);
    ylabel(hca,{'L_{2}/L_{3}'})
    %irf_zoom(hca,'y',[0 1])
end
if 0 % <|B|> 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,mB?);',sc);
    ylabel(hca,{'<|B|>','[nT]'})
end
if 0 % E/B ratio
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,irf_tappl(EB,'*1e-3'));
    c_eval('irf_plot(hca,irf_tappl(EB?,''*1e-3''));',sc);
    ylabel(hca,{'<|E|>/<|B|>','[km/s]'})
    %irf_zoom(hca,'y',[0 50000])
end

irf_zoom(h,'x',tint_a+[+3 0])
irf_plot_axis_align

%% make quiver plot of dE and dB
tint_q = toepoch([2002 03 30 13 11 45.80; 2002 03 30 13 11 46.50])';
sc = 2;
ffilt=10;
c_eval('Eq = irf_tlim(facE?,tint_q);',sc);
c_eval('Bq = irf_filt(irf_tlim(facB?(:,1:4),tint_q),10,0);',sc);

nq = size(Eq,1);
topl = 1:2:nq;

v=400; % km/s
T = Eq(topl,1);
X = (T-T(1))*v;
Y = X*0;


npl=3;
for kk=1:npl; h(kk)=subplot(npl,1,kk); end
isub=1;
if 1
    hca=h(isub);isub=isub+1;
    hold(hca,'on');
    quiver(hca,X,Y,detrend(Eq(topl,4)),detrend(Eq(topl,3)),'color',[0.2 0.8 0.8]); 
    quiver(hca,X,Y,detrend(Bq(topl,4)),detrend(Bq(topl,3)),'color',[0.8 0.2 0.8])
    %irf_legend({'E','B'},[1 0.95])
    legend(hca,'E','B')    
    bpar= find(detrend(Bq(topl,2))>0);
    bapar= find(detrend(Bq(topl,2))<0);
    plot(hca,X(bpar),Y(bpar),'o',X(bapar),Y(bapar),'x','markersize',10,'color',[0.8 0.2 0.8])
    set(hca,'ylim',20*[-1 1])
    %axis(hca,'square')
    
    %axis(hca,'equal')
end
if 1
    hca=h(isub);isub=isub+1;
    %legend(hca,'E','B')
    irf_plot(hca,Eq)
    irf_legend(hca,{'L','M','N'},[1 0.95])
    irf_zoom(hca,'x',tint_q)
end
if 1
    hca=h(isub);isub=isub+1;
    %legend(hca,'E','B')
    irf_plot(hca,[Bq(:,1) detrend(Bq(:,2)) detrend(Bq(:,3)) detrend(Bq(:,4))])
    irf_legend(hca,{'L','M','N'},[1 0.95])
    irf_zoom(hca,'x',tint_q)
end

%% compare the EB polarizations in a local fiald aligned MVA system
tint_q = toepoch([2002 03 30 13 11 45.80; 2002 03 30 13 11 46.50])';
tint_q = toepoch([2002 03 30 13 11 49.65; 2002 03 30 13 11 49.85])';
sc = 2;
ffilt=10;
c_eval('Eq = irf_tlim(facE?,tint_q);',sc);
c_eval('Bq = irf_filt(irf_tlim(facB?(:,1:4),tint_q),20,0);',sc);
[~,l,v] = irf_minvar(Bq);
[~,lE,vE] = irf_minvar(Eq);

c_eval('d1 = cn_average2(irf_tlim(facB?(:,1:4),tint_q),1)/norm(cn_average2(irf_tlim(facB?(:,1:4),tint_q),1));',sc); % B
d3 = cross(d1,cross(d1,v(3,:))); % v3, but perp to B
d2 = cross(d3,d1); % third direction
fBq = irf_lmn(Bq,d1,d2,d3);
fEq = irf_lmn(Eq,d1,d2,d3);

%
h=irf_plot(3);
isub=1;
coord_leg = {'L','M','N'};
sc = 2;
if 1 % background magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_abs(facB?(:,1:4)))',sc);
    ylabel(hca,'B [nT]')
    irf_legend(hca,{'||','N','MVA'},[1 0.95])
end
if 1 % MVA electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,fEq)',sc);
    ylabel(hca,'E [mV/m]') 
    irf_legend(hca,{'||','N','MVA'},[1 0.95]) 
end
if 1 % MVA wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,fBq)',sc);
    ylabel(hca,'\delta B [nT]') 
    irf_legend(hca,{'||','N','MVA'},[1 0.95]) 
end
if 0 % FAC electric field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,facE?)',sc);
    ylabel(hca,'E_{fac} [mV/m]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end
if 0 % FAC wave magnetic field
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_filt(facB?(:,1:4),30,0))',sc);
    ylabel(hca,'\delta B_{fac} [nT]') 
    irf_legend(hca,{'L','M','N'},[1 0.95]) 
end
if 0 % k directions from bool and tool, 3 panels
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 2]),kmvas?(:,[1 2])},''comp'')',sc);
    ylabel(hca,'k_L') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 3]),kmvas?(:,[1 3])},''comp'')',sc);
    ylabel(hca,'k_M') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,{ks?(:,[1 4]),kmvas?(:,[1 4])},''comp'')',sc);
    ylabel(hca,'k_N') 
    set(hca,'ylim',1.1*[-1 1])
    irf_legend(hca,{'tool','bool'},[1 0.95]) 
end
if 0 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l12_?);',sc);
    ylabel(hca,{'L_{1}/L_{2}'})
    %irf_zoom(hca,'y',[0 1])
end
if 0 % eigenvalue ratios
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,l23_?);',sc);
    ylabel(hca,{'L_{2}/L_{3}'})
    %irf_zoom(hca,'y',[0 1])
end

irf_zoom(h,'x',tint_q)
irf_plot_axis_align