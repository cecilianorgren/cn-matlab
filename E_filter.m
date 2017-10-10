tic
figure(6)
h=irf_plot(12);

%% Non-filtered
dt=diE3zoom(2,1)-diE3zoom(1,1);
fs=1/dt;
overlap=0;
nfft=256;

c_eval('E?fft=irf_powerfft(diE?zoom,nfft,fs,overlap);',3:4);
irf_plot(h(1),diE3zoom(:,1:3));
irf_spectrogram([h(2) h(3)],E3fft);
irf_plot(h(4),diE4zoom(:,1:3));
irf_spectrogram([h(5) h(6)],E4fft);

%% Filtered
if 0
E_filt=caa_filter_e(diE3zoom(:,1:3),0.05); % 20Hz?
irf_plot(h(4),E_filt);
Efft=irf_powerfft(E_filt,nfft,fs,overlap);
irf_spectrogram([h(5) h(6)],Efft);
end
%% Other filter

c_eval('E?_filt=irf_filt(diE?zoom(:,1:3),10,180,fs,3);',3:4);
c_eval('E?fft=irf_powerfft(E?_filt,nfft,fs,overlap);',3:4);

irf_plot(h(7),E3_filt);
irf_spectrogram([h(8) h(9)],E3fft);

irf_plot(h(10),E4_filt);
irf_spectrogram([h(11) h(12)],E4fft);



%%
irf_zoom(h,'x',tint_zoom(5,:));
 irf_plot_axis_align
    add_timeaxis(h(6),'usefig');
    
 set(gcf,'PaperPositionMode','auto')   
 print -dpdf filter.pdf
    toc