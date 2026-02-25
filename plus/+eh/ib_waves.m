% Looks at wave spectra during the two internal burst intervals available.
cd /Users/Cecilia/Data/BM/20070831;
sclist=3:4;
c_eval('EB? = c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sclist);
load matlabB
%% make plasma parameters
c_eval('tint?=[EB?(1,1) EB?(end,1)];',sclist);
tint=[min([EB3(1,1) EB4(1,1)]) max([EB3(end,1) EB4(end,1)])];
%c_eval('B?=irf_tlim(diB?,tint);',sclist); % zoomed B
%B=mean([B3(:,5);B4(:,5)]);
B=irf_add(0.5,irf_tlim(diB3,tint),0.5,irf_tlim(diB4,tint));

fce=irf_plasma_calc(B,0.1,0,100,100,'Fce');
fpe=irf_plasma_calc(B,0.1,0,100,100,'Fpe');
fcp=irf_plasma_calc(B,0.1,0,100,100,'Fcp');
fpp=irf_plasma_calc(B,0.1,0,100,100,'Fpp');
flh=irf_plasma_calc(B,0.1,0,100,100,'Flh');
%% making power spectrum
c_eval('notNan?=~isnan(EB?(:,2));',sclist);
nfft=256; fs=1/diff(EB3(1:2,1)); overlap=50;
c_eval('fftEB?_256 = irf_powerfft(EB?(notNan?,:),nfft,fs);',sclist)
c_eval('wletEB?_100 = irf_wavelet(EB?(notNan?,:),''nf'',100,''fs'',fs);',sclist)
nfft=512; fs=1/diff(EB3(1:2,1)); overlap=50;
c_eval('fftEB?_512 = irf_powerfft(EB?(notNan?,:),nfft,fs);',sclist)
nfft=1024; fs=1/diff(EB3(1:2,1)); overlap=50;
c_eval('fftEB?_1024 = irf_powerfft(EB?(notNan?,:),nfft,fs);',sclist)
%% plot
h=irf_plot(6);
iSub=1;
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB3_256); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    text(hca,1,fce(1,2)*1e-3,'f_{ce}');
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB4_256); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB3_512); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB4_512); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB3_1024); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,fftEB4_1024); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end
irf_zoom(h,'x',tint)

%% more info plots
fftEB3.plot_type='log';
h=irf_plot(3);
iSub=1;
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_spectrogram(hca,fftEB3_256); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    %irf_plot(hca,irf_tappl(fcp,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpp,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(flh,'*1e-3')); % kHz
    %text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log','ylim',[1e-2 5e0])
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_spectrogram(hca,wletEB3_100); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
    %irf_plot(hca,irf_tappl(fcp,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpp,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(flh,'*1e-3')); % kHz
    %text(hca,1,fce(1,2)*1e-3,'f_{ce}');
    set(hca,'yscale','log','ylim',[1e-2 5e0])
end
if 1
    hca=h(iSub); iSub=iSub+1;
    irf_plot(hca,EB3); hold(hca,'on');
    irf_plot(hca,irf_tappl(fce,'*1e-3')); % kHz
    irf_plot(hca,irf_tappl(fpe,'*1e-3')); % kHz
end


irf_plot_axis_align
irf_zoom(h,'x',tint3)