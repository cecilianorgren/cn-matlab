%%
tint_zoom=tint(9,:);

%%
tint_zoom=[toepoch([2007 08 31 10 15 00]) toepoch([2007 08 31 10 25 00])];
dt=diE3zoom(2,1)-diE3zoom(1,1);
fs=1/dt;
overlap=0;
nfft=256;
i_start=find(diE3(:,1)>tint_zoom(1),1);
i_end=find(diE3(:,1)>tint_zoom(2),1);
diE3zoom=diE3(i_start-nfft:i_end+nfft,:);

diE3fft=irf_powerfft(diE3zoom,nfft,fs,overlap);
%%
sizeNi3=size(Ni3,1);
Bindex=fix(linspace(1,size(absB3,1),sizeNi3));
absB3red=absB3(Bindex,[1 5]);
%%
%index=0;
for t=1:sizeNi3
    %index=index+1;
    %disp([num2str(index)]);
    flh(t,1)=absB3red(t,1);
    flh(t,2)=irf_plasma_calc(absB3red(t,2),Ni3(t,2),0,2500,Ti3eV(t),'Flh');
end
%%

h=irf_plot(9);

irf_plot(h(1),absB3(:,[1 5]));
irf_plot(h(2),absB3red);
irf_plot(h(3),diE3);
irf_plot(h(4),Ni3);
irf_plot(h(5),Ti3eV);
irf_plot(h(6),Te3pereV);
irf_plot(h(7),flh);
irf_spectrogram([h(8) h(9)],diE3fft); hold(h(8),'on'); hold(h(9),'on');
        ylabel(h(8),'f [Hz]');
        ylabel(h(9),'f [Hz]');
        irf_legend(h(8),'C3_X',[0.02 0.90]);
        irf_legend(h(9),'C3_Y',[0.02 0.90]);

irf_plot(h(8),flh);
set(h(8),'yscale','log');
set(h(8),'ytick',[1 1e1 1e2]);
irf_plot(h(9),flh);
set(h(9),'yscale','log');
set(h(9),'ytick',[1 1e1 1e2]);

ylabel(h(1),'|B| [nT]')
ylabel(h(2),'|B|_{red} [nT]')
ylabel(h(3),'E [mV/m]')
ylabel(h(4),'N_i [cm^{-3}]')
ylabel(h(5),'T_i [nT]')
ylabel(h(6),'T_{e,\perp} [nT]')
ylabel(h(7),'f_{lh} [nT]')
ylabel(h(8),'f [Hz]');
ylabel(h(9),'f [Hz]');

irf_zoom(h,'x',tint_zoom)
irf_plot_axis_align
set(gcf,'PaperPositionMode','auto','Position',position); 
print -dpdf 20070831_spectrum_flh_overlap.pdf
