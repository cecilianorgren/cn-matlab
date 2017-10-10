%% Ion data

c_eval('ipaddist? = mms.get_pitchangledist(disDist?,dmpaB?brst);')

%%
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

energyLevels = [5 10 15 20 25 30];
ic=1:4;
for ii = 1:numel(energyLevels)
c_eval('ipad_spec?_!.p = squeeze(ipaddist?.data(:,energyLevels(ii),:)); ipad_spec?_!.t = ipaddist?.time.epochUnix; ipad_spec?_!.f = 15:15:180; ipad_spec?_!.energy = energy(!);',ic,energyLevels(ii))
end
%% Overview
tint = irf.tint('2015-10-16T10:33:00.00Z/2015-10-16T10:34:20.00Z'); % magnetosphere-magnetosheath-magnetosphere
h = irf_plot(12);

ic  =4;
% Magnetic field
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.7]);



hca = irf_panel('brst n');
set(hca,'ColorOrder',mms_colors('xyza'))
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'n_e','n_i'},[1.01 0.5]);

hca = irf_panel('brst vi');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{vi?_lowres.tlim(tint).x,vi?_lowres.tlim(tint).y,vi?_lowres.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'v_i','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[1.01 0.7]);

hca=irf_panel('idist');
c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{i,OMNI}','(eV)'},'Interpreter','tex');

for ii = 1:numel(energyLevels)
  hca=irf_panel(irf_ssub('idist_pad ?',energyLevels(ii)));  
  c_eval('irf_spectrogram(hca,ipad_spec?_!,''log'',''donotfitcolorbarlabel'');',ic,energyLevels(ii))
  %set(hca,'yscale','log');
  set(hca,'ytick',15:30:180);
  ylabel(hca,{'\theta','(?)'},'Interpreter','tex');
  irf_legend(hca,['E_i = ' num2str(energy(energyLevels(ii)),'%.0f') ' eV',],[0.02 0.98])
end
if 0
  hca=irf_panel('idist_pad 2');
  c_eval('irf_spectrogram(hca,ipad_spec?_20,''log'',''donotfitcolorbarlabel'');',ic)
  %set(hca,'yscale','log');
  set(hca,'ytick',15:30:180);
  ylabel(hca,{'\theta','(?)'},'Interpreter','tex');
end


hca = irf_panel('j curl');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{jbrst.x,jbrst.y,jbrst.z},'comp');
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'J_x','J_y','J_z'},[1.01 0.7]);

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{e,OMNI}','(eV)'},'Interpreter','tex');


irf_zoom(h,'x',irf.tint('2015-10-16T10:33:00.00Z/2015-10-16T10:34:20.00Z'))
%irf_zoom(h([1:5 10:11]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
%labelling
labels = {'a','b','c','d','e','f','g','h','j','k','l','m','n','o'};
for ii = 1:numel(h);
  irf_legend(h(ii),labels{ii},[0.98 0.98],'color',[0 0 0])
end
