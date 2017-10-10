%% Load and display/plot whisper data
cd /Users/Cecilia/Data/BM/20070831;
%% C*_CP_WHI_WAVE_FORM_ENERGY
[varEWFPD,dobjEWFPD,matEWFPD,unitEWFPD]=c_caa_var_get('Electric_Wave_Form_Power_Density__C4_CP_WHI_WAVE_FORM_ENERGY');
figure(1);
irf_plot('Electric_Spectral_Power_Density__C4_CP_WHI_PASSIVE_ACTIVE')

%%

%% C*_CP_WHI_PASSIVE_ACTIVE
[varESPDap,dobjESPDap,matESPDap,unitESPDap]=c_caa_var_get('Electric_Spectral_Power_Density__C4_CP_WHI_PASSIVE_ACTIVE');
%figure(1);
%irf_plot('Electric_Spectral_Power_Density__C4_CP_WHI_PASSIVE_ACTIVE')

%%
[varESPD,dobjESPD,matESPD,unitESPD]=c_caa_var_get('Electric_Spectral_Power_Density__C4_CP_WHI_ACTIVE');

%%
[varNAT,dobjNAT,matNAT,unitNAT]=c_caa_var_get('Electric_Spectral_Power_Density__C4_CP_WHI_NATURAL');
freqNAT=c_caa_var_get('Spectral_Frequencies__C4_CP_WHI_NATURAL','mat');
freqUnits=c_caa_var_get('Spectral_Frequencies__C4_CP_WHI_NATURAL','units');
tint = toepoch([2007 08 31 10 00 00;2007 08 31 10 25 00])';
[limData,ind] = irf_tlim(matNAT(:,1),tint);
units=irf_units;
fpe = [peaNe4(:,1) sqrt(units.e^2*peaNe4(:,2)*1e6/units.me/units.eps0)/1000/2/pi];
fce = [diB4(:,1) units.e*diB4(:,5)*1e-9/units.me/1000/2/pi];
fuh = irf_add(1,fpe,1,fce);
%fpe2 = [peaNe4(:,1) 9e3*sqrt(peaNe4(:,2))/1000];

specrec.t = matNAT(ind,1);
specrec.f = freqNAT;
specrec.p = matNAT(ind,[2:end]);
specrec.f_label = 'f [kHz]';
specrec.p_label = unitNAT;
specrec.plot_type = 'log';

h = irf_spectrogram(specrec);   
set(h,'yscale','log')
title(h,'Electric Spectral Power Density')
hold(h,'on');
irf_plot(h,fpe,'k')
irf_plot(h,fce,'k')
irf_plot(h,fuh,'k')
text(80,18,'f_{pe}','fontsize',14)
text(14*60,4,'f_{ce}','fontsize',14)
text(14*60,15,'f_{uh}','fontsize',14)
hold(h,'off')


