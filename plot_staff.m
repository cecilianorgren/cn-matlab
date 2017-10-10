%% Load and display/plot whisper data
cd /Users/Cecilia/Data/BM/20070831;
if ~exist('peaNe4','var'), load matlabN; end
if ~exist('diB4','var'), load matlabdiEB; end
%% C*_CP_STA_PPP
[varSTAPSD,dobjSTAPSD,matSTAPSD,unitSTAPSD]=c_caa_var_get('EE_xxyy_sr2__C4_CP_STA_PSD');
freqSTAPSD = c_caa_var_get('Frequency__C4_CP_STA_PSD','mat');
freqUnits = c_caa_var_get('Frequency__C4_CP_STA_PSD','units');
tint = toepoch([2007 08 31 10 00 00;2007 08 31 10 25 00])';
[limData,ind] = irf_tlim(matSTAPSD(:,1),tint);
units=irf_units;
fpe = [peaNe4(:,1) sqrt(units.e^2*peaNe4(:,2)*1e6/units.me/units.eps0)/2/pi];
fce = [diB4(:,1) units.e*diB4(:,5)*1e-9/units.me/2/pi];
fuh = irf_add(1,fpe,1,fce);

specrecx.t = matSTAPSD.t(ind,1);
specrecx.f = freqSTAPSD;
specrecx.p = matSTAPSD.data(ind,:,1); % xx-component
specrecx.f_label = ['f [' freqUnits ']'];
specrecx.p_label = unitSTAPSD;
specrecx.plot_type = 'log';

specrecy.t = matSTAPSD.t(ind,1);
specrecy.f = freqSTAPSD;
specrecy.p = matSTAPSD.data(ind,:,2); % xx-component
specrecy.f_label = ['f [' freqUnits ']'];
specrecy.p_label = unitSTAPSD;
specrecy.plot_type = 'log';

structs={specrecx,specrecy};
%
h = irf_plot(3,'newfigure',77);

for k=1:2;
irf_spectrogram(h(k),structs{k});   
set(h(k),'yscale','log')
hold(h(k),'on');
irf_plot(h(k),fpe,'k')
irf_plot(h(k),fce,'k')
irf_plot(h(k),fuh,'k')
text(12.5*60,1.5e3,'f_{pe}','fontsize',14)
text(14*60,4.5e2,'f_{ce}','fontsize',14)
text(14*60,2e3,'f_{uh}','fontsize',14)
hold(h(k),'off')
end
title(h(1),varSTAPSD.FIELDNAM)

if 1
    hca=h(3);
    irf_plot(hca,diE4)
    ylabel('E4 [mV/m] ISR2')    
end
irf_plot_axis_align
irf_zoom(h,'x',toepoch([2007 08 31 10 12 00;2007 08 31 10 25 00])')

%% C*_CP_STA_PPP : magnetic field
[varSTAPSD,dobjSTAPSD,matSTAPSD,unitSTAPSD]=c_caa_var_get('BB_xxyyzz_sr2__C4_CP_STA_PSD');
freqSTAPSD = c_caa_var_get('Frequency__C4_CP_STA_PSD','mat');
freqUnits = c_caa_var_get('Frequency__C4_CP_STA_PSD','units');
tint = toepoch([2007 08 31 10 00 00;2007 08 31 10 25 00])';
[limData,ind] = irf_tlim(matSTAPSD(:,1),tint);
units=irf_units;
fpe = [peaNe4(:,1) sqrt(units.e^2*peaNe4(:,2)*1e6/units.me/units.eps0)/2/pi];
fce = [diB4(:,1) units.e*diB4(:,5)*1e-9/units.me/2/pi];
fuh = irf_add(1,fpe,1,fce);

specrecx.t = matSTAPSD.t(ind,1);
specrecx.f = freqSTAPSD;
specrecx.p = matSTAPSD.data(ind,:,1); % xx-component
specrecx.f_label = ['f [' freqUnits ']'];
specrecx.p_label = unitSTAPSD;
specrecx.plot_type = 'log';

specrecy.t = matSTAPSD.t(ind,1);
specrecy.f = freqSTAPSD;
specrecy.p = matSTAPSD.data(ind,:,2); % yy-component
specrecy.f_label = ['f [' freqUnits ']'];
specrecy.p_label = unitSTAPSD;
specrecy.plot_type = 'log';

specrecz.t = matSTAPSD.t(ind,1);
specrecz.f = freqSTAPSD;
specrecz.p = matSTAPSD.data(ind,:,3); % zz-component
specrecz.f_label = ['f [' freqUnits ']'];
specrecz.p_label = unitSTAPSD;
specrecz.plot_type = 'log';

structs={specrecx,specrecy,specrecz};
%
h = irf_plot(4,'newfigure',77);

for k=1:3;
irf_spectrogram(h(k),structs{k});   
set(h(k),'yscale','log')
set(h(k),'ylim',[50 4500])
hold(h(k),'on');
irf_plot(h(k),fpe,'k')
irf_plot(h(k),fce,'k')
irf_plot(h(k),fuh,'k')
text(12.5*60,1.5e3,'f_{pe}','fontsize',14)
text(14*60,4.5e2,'f_{ce}','fontsize',14)
text(14*60,2e3,'f_{uh}','fontsize',14)
hold(h(k),'off')
end
clim=get(h(1),'clim')
set(h(2),'clim',clim)
set(h(3),'clim',clim)
for k=1:3; set(h(k),'clim',[-13 -5.5]); end

title(h(1),varSTAPSD.FIELDNAM)

if 1
    hca=h(4);
    irf_plot(hca,diE4)
    ylabel('E4 [mV/m] ISR2')    
end
irf_plot_axis_align
irf_zoom(h,'x',toepoch([2007 08 31 10 12 00;2007 08 31 10 25 00])')