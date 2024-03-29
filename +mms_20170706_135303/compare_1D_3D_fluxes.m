units = irf_units;
tint_zoom = irf.tint('2017-07-06T13:54:00.000Z/2017-07-06T13:54:10.000Z');
tint_zoom = irf.tint('2017-07-06T13:53:55.000Z/2017-07-06T13:54:15.000Z');

% EDI energy and corresponding velocity
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

v_edi_edges = [v_edi_minus,v_edi_plus]*1e-3;
pa_edi_edges = -45:11.25:45;
pa_edi_centers = [(-45+11.25):11.25:0 0:11.25:(45-11.25)];
az_edges = 0:10:360;

tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
tint_zoom = tint_zoom + [-1 1];

c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',1:4)
c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4);
c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
c_eval('ePitch?_flux_fpi_bin_closest_to_edi = ePitch?_flux_fpi.elim(500);',1:4);
%c_eval('ePitch?_flux_fpi_edi_range = ePitch?_flux_fpi.rebin({});',1:4);

%c_eval('eRed?_fpi_edi_range = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',(-[v_edi_minus v_edi_plus])*1e-3);',1:4)
c_eval('eRed?_fpi_edi_range_180 = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg_edges'',([-v_edi_plus -v_edi_minus ])*1e-3);',1)
c_eval('eRed?_fpi_flux_edi_range = eRed?_fpi_edi_range_180.flux_red;',1)

%% Plot, compare EDI and FPI for different spacecraft
npanels = 5;
h = irf_plot(npanels);
isub = 1;
elim = 500;

hlink = gobjects(0);
% for sc = 1:4 % EDI vs FPI
%   c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',sc)
%   c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',sc)  
% end

if 1 % Epar
  hca = irf_panel('Epar');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp')
  irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end

palim = 180-12*3;
for ic = 1:4
if 1 % j FPI EDI 180 
  hca = irf_panel(sprintf('mms %g, j_%g',ic,palim));      
  c_eval('pa_edi = ePitch?_flux_edi.palim(palim).depend{2};',ic)
  c_eval('en_edi = ePitch?_flux_edi.palim(palim).depend{1}(1);',ic)
  c_eval('pa_fpi = ePitch?_flux_fpi.elim(elim).palim(palim).depend{2};',ic)
  c_eval('en_fpi = ePitch?_flux_fpi.elim(elim).palim(palim).depend{1}(1);',ic)  
  
  c_eval('j_edi = ePitch?_flux_edi.palim(palim)*1e-6;',ic)
  c_eval('j_fpi = ePitch?_flux_fpi.elim(elim).palim(palim)*1e-6;',ic)  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{j_edi,j_fpi},'comp');  
  irf_legend(hca,{'EDI','FPI'},[0.98 0.98])
  hca.YLabel.String = {sprintf('j_{MMS %g}',ic);'(10^{6} cm^{-2} s^{-1} sr^{-1})'};
  irf_legend(hca,{['E_{FPI/EDI} = ' sprintf('%.2f',en_fpi) '/' sprintf('%.2f',en_edi)];...
    ['\theta_{FPI/EDI} = ' sprintf('%.2f',pa_fpi) '/' sprintf('%.2f',pa_edi)]},[0.02 0.98],'color',[0 0 0])
  
  hl = findall(hca.Children,'type','line');
  hl(1).LineWidth = 1;
  
end
end

%irf_zoom(h,'x',tint_zoom+[-1 1])
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
set(irf_panel('Epar'),'ylim',[-60 60]);
c_eval('h(?).YLim = [0 10];',2:5)

%cn.print(sprintf('j_EDI_FPI_pa_%.0f_ylim10_zoomout',pa_edi))
cn.print(sprintf('j_EDI_FPI_pa_%.0f_ylim10_zoomout_all',pa_edi))


%% Plot compare flux at different pitch angles
npanels = 8;
h = irf_plot(npanels);
elim = 500;
eflow = 10; efhigh = 0;
timeline = gseVe1.time;

if 0 % B
  hca = irf_panel('B1');
  irf_plot(hca,{gseB1})
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hca.YLabel.String = {'B_1','(nT)'};
end
if 1 % Epar
  hca = irf_panel('Epar');  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp')
  irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{||}','(mV/m)'};
end
if 1 % Eperpx
  hca = irf_panel('Eperpx');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.x,gseE2perp.x,gseE3perp.x,gseE4perp.x},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{\perp,x}','(mV/m)'};
end
if 1 % Eperpy
  hca = irf_panel('Eperpy');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.y,gseE2perp.y,gseE3perp.y,gseE4perp.y},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};
end
if 0 % Eperpy, filt
  hca = irf_panel('Eperpy');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.y.filt(eflow,efhigh),gseE2perp.y.filt(eflow,efhigh),gseE3perp.y.filt(eflow,efhigh),gseE4perp.y.filt(eflow,efhigh)},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};
end
if 1 % Eperpz
  hca = irf_panel('Eperpz');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1perp.z,gseE2perp.z,gseE3perp.z,gseE4perp.z},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'E_{\perp,z}','(mV/m)'};
end
if 0 % EDI 180 
  hca = irf_panel('j EDI 180-0*12');    
  palim = 180-12*0;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,...
                ePitch2_flux_edi.palim(palim)*1e-6,...
                ePitch3_flux_edi.palim(palim)*1e-6,...
                ePitch4_flux_edi.palim(palim)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 0 % EDI 180 
  hca = irf_panel('j EDI 180-1*12');    
  palim = 180-12*1;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,...
                ePitch2_flux_edi.palim(palim)*1e-6,...
                ePitch3_flux_edi.palim(palim)*1e-6,...
                ePitch4_flux_edi.palim(palim)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 0 % EDI 180 
  hca = irf_panel('j EDI 180-2*12');    
  palim = 180-12*2;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,...
                ePitch2_flux_edi.palim(palim)*1e-6,...
                ePitch3_flux_edi.palim(palim)*1e-6,...
                ePitch4_flux_edi.palim(palim)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 0 % EDI 180 
  hca = irf_panel('j EDI 180-3*12');    
  palim = 180-12*3;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,...
                ePitch2_flux_edi.palim(palim)*1e-6,...
                ePitch3_flux_edi.palim(palim)*1e-6,...
                ePitch4_flux_edi.palim(palim)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 1 % EDI 180 resample fpi
  hca = irf_panel('j EDI 180-0*12');    
  palim = 180-12*0;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch2_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch3_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch4_flux_edi.palim(palim).resample(timeline)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 1 % EDI 180 
  hca = irf_panel('j EDI 180-1*12');    
  palim = 180-12*1;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch2_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch3_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch4_flux_edi.palim(palim).resample(timeline)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 1 % EDI 180 
  hca = irf_panel('j EDI 180-2*12');    
  palim = 180-12*2;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch2_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch3_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch4_flux_edi.palim(palim).resample(timeline)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
if 1 % EDI 180 
  hca = irf_panel('j EDI 180-3*12');    
  palim = 180-12*3;
  c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',1)
  c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',1)  
  
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch2_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch3_flux_edi.palim(palim).resample(timeline)*1e-6,...
                ePitch4_flux_edi.palim(palim).resample(timeline)*1e-6},'comp')
   irf_legend(hca,{'MMS 1','MMS 2','MMS 3','MMS 4'}',[1.01 0.98])
  hca.YLabel.String = {'j^{EDI}';'(10^{6} cm^{-2} s^{-1})'};
  irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
end
%irf_zoom(h,'x',tint_zoom+[-1 1])
irf_zoom(h,'y')

%c_eval('h(?).YLim = [0 8];',5:8)
c_eval('h(?).YLim = [0 9];',5:8)
%% Plot, compare fpi to edi for individual spacecraft
for sc = 1:4
  for palim = 180-[0 1 2 3]*12
    npanels = 4;
    h = irf_plot(npanels);
    isub = 1;
    %sc = 4;
    legends_edi = {};
    for ipa=1:8 legends_edi{ipa} = sprintf('[%.2f %.2f]',ePitch1_flux_edi.depend{2}(ipa)-ePitch1_flux_edi.ancillary.delta_pitchangle_minus(ipa)*0.5,ePitch1_flux_edi.depend{2}(ipa)+ePitch1_flux_edi.ancillary.delta_pitchangle_plus(ipa)*0.5); end

    %palim = 180-12;
    elim = 500;
    c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',sc)
    c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',sc)


    if 0 % edi flux par
      hca = h(isub); isub = isub + 1;
      palim = [0 45];
      irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6})
      irf_legend(hca,legends_edi(1:4),[0.01 0.98])
      hca.YLabel.String = {'Flux edi','(10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 0 % edi fpi 0
      hca = h(isub); isub = isub + 1;
      palim = 0;
      irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6},'comp')
      irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
      hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 0 % edi flux apar
      hca = h(isub); isub = isub + 1;
      palim = [135 180];
      irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6})
      irf_legend(hca,legends_edi(5:8),[0.01 0.98])
      hca.YLabel.String = {'Flux edi','(10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 1 % edi fpi 180
      hca = h(isub); isub = isub + 1;
      set(hca,'ColorOrder',pic_colors('xy'))
      c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6,ePitch?_flux_fpi.palim(palim).elim(elim)*1e-6},''comp'');',sc)
      irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
      hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 1 % edi fpi 180, edi resampled to fpi
      hca = h(isub); isub = isub + 1;
      set(hca,'ColorOrder',pic_colors('xy'))
      c_eval('irf_plot(hca,{ePitch?_flux_edi.resample(ePitch1_flux_fpi).palim(palim)*1e-6,ePitch?_flux_fpi.palim(palim).elim(elim)*1e-6},''comp'');',sc)
      irf_legend(hca,{'EDI','FPI'},[0.01 0.98])
      hca.YLabel.String = {'Flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 1 % fpi red edi e range
      hca = h(isub); isub = isub + 1;
      set(hca,'ColorOrder',pic_colors('z'))
      c_eval('irf_plot(hca,{eRed?_fpi_edi_range_180.flux_red*1e-6},''comp'');',sc)
      c_eval('irf_legend(hca,{sprintf(''v = [%.0f, %.0f] km/s'',eRed?_fpi_edi_range_180.ancillary.v_edges(1)*1e-3,eRed?_fpi_edi_range_180.ancillary.v_edges(2)*1e-3)},[0.01 0.05]);',sc)
      hca.YLabel.String = {'Reduced flux','FPI','(10^6 s^{-1}cm^{-2})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 1 % fpi red edi e range, edi fpi 180
      hca = h(isub); isub = isub + 1;      
      set(hca,'ColorOrder',pic_colors('xyz'))      
      c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6,ePitch?_flux_fpi.palim(palim).elim(elim)*1e-6,eRed?_fpi_edi_range_180.flux_red*1e-6},''comp'');',sc)
      c_eval('v1 = eRed?_fpi_edi_range_180.ancillary.v_edges(1)*1e-3;',sc)
      c_eval('v2 = eRed?_fpi_edi_range_180.ancillary.v_edges(2)*1e-3;',sc)
      irf_legend(hca,{'EDI fov',...
        'FPI fov',...
        sprintf('FPI red: v = [%.0f, %.0f] km/s',v1,v2)},...
        [0.01 0.98])
      hca.YLabel.String = {'Flux','red (10^6 s^{-1}cm^{-2})','fov (10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    if 0 % fpi red edi e range apar, edi fpi 0, just compare for debugging
      hca = h(isub); isub = isub + 1;
      palim = 0;
      set(hca,'ColorOrder',pic_colors('xyz'))
      irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch1_flux_fpi.palim(palim).elim(500)*1e-6,eRed1_fpi_edi_range_180.flux_red*1e-6},'comp')
      irf_legend(hca,{'EDI fov',...
        'FPI fov',...
        sprintf('FPI red: v = [%.0f, %.0f] km/s',eRed1_fpi_edi_range_180.ancillary.v_edges(1)*1e-3,eRed1_fpi_edi_range_180.ancillary.v_edges(2)*1e-3)},...
        [0.01 0.98])
      hca.YLabel.String = {'Flux','red (10^6 s^{-1}cm^{-2})','fov (10^6 s^{-1}cm^{-2}sr^{-1})'};
      hca.YLabel.Interpreter = 'tex';
    end
    h(1).Title.String = sprintf('MMS %g, pitch angle = %.2f, energy = %.2f',sc,pa,en);

    irf_zoom(h,'x',tint)
    drawnow
    %cn.print(sprintf('edi_fpi_comp_mms%g_en=%.0f_pa=%.0f_zoom_0',sc,en,pa))
    
    irf_zoom(h,'x',tint_zoom+[-5 5])
    drawnow
    %cn.print(sprintf('edi_fpi_comp_mms%g_en=%.0f_pa=%.0f_zoom_1',sc,en,pa))
    
    irf_zoom(h,'x',tint_zoom)
    drawnow
    %cn.print(sprintf('edi_fpi_comp_mms%g_en=%.0f_pa=%.0f_zoom_2',sc,en,pa))
  end
end

%% Plot. bincounts(edi1,fpi1), bincounts(edi2,fpi2)
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
palim = 180-12*0;
elim = 500;

hlink = gobjects(0);
for sc = 1:4 % EDI vs FPI
  if 1 % EDI vs FPI 
    hca = h(isub); isub = isub + 1;    
    hlink(end+1) = hca;
    c_eval('pa = ePitch?_flux_edi.palim(palim).depend{2};',sc)
    c_eval('en = ePitch?_flux_edi.palim(palim).depend{1}(1);',sc)
    if 1 % downsample EDI to FPI
      c_eval('xx = ePitch?_flux_edi.palim(palim).resample(ePitch?_flux_fpi.time).data*1e-6;',sc)
      c_eval('yy = ePitch?_flux_fpi.palim(palim).elim(elim).data*1e-6;',sc)
      clim = [0 2];
    else % upsample FPI to EDI
      c_eval('xx = ePitch?_flux_edi.palim(palim).data*1e-6;',sc)
      c_eval('yy = ePitch?_flux_fpi.palim(palim).resample(ePitch?_flux_edi.time).elim(elim).data*1e-6;',sc)
      clim = [0 3];
    end
    edges = linspace(0,10,100);
    data = [xx yy];
    [N edges mid loc] = histcn(data,edges,edges);
    N(N==0) = NaN;
    pcolor(hca,mid{1:2},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10} counts';
    colormap(hca,pic_colors('pasteljet'))
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    axis(hca,'square')
    hca.CLim = clim;
    %plot(hca,xx,yy,'.')
    hca.XTick = hca.YTick;
    hca.XLabel.String = sprintf('EDI (%s 10^6)',ePitch1_flux_edi.units);
    hca.YLabel.String = sprintf('FPI (%s 10^6)',ePitch1_flux_fpi.units);
    hca.Title.String = sprintf('MMS %g',sc);
    irf_legend(hca,{['energy = ' sprintf('%.2f',en)];['\theta = ' sprintf('%.2f',pa)]},[0.02 0.98],'color',[0 0 0])
    hold(hca,'on')
    plot(hca,hca.XLim,hca.XLim,'color',[0 0 0])
    hold(hca,'off')
  end  
end
hlinks = linkprop(hlink,{'XLim','YLim','CLim'});
hca.YLim = [0 10];
hca.XLim = [0 10];
c_eval('h(?).XTick = h(?).YTick;',1:4)


%%
%c_eval('eDred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)

%c_eval('eFlred? = ePDist?.tlim(tint_zoom).reduce(''1D'',dmpaB?.resample(ePDist?.tlim(tint_zoom)).norm,''vg'',([v_edi_minus v_edi_plus]+0.25*v_edi_plusminus*[1 -1])*1e-3);',1:4)
%dv = v_edi_plusminus; % m/s
%c_eval('eDred?_ = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)); eDred?_.units = eDred?.units;',1:4)
%c_eval('eFlux_red? = irf.ts_scalar(eDred?.time,sum(eDred?.data,2)*v_edi*dv); eFlux_red?.units = ''1/s/m^2'';',1:4)
%c_eval('eFlux_red? = eFlux_red?*1e-4; eFlux_red?.units = ''1/s/cm^2'';',1:4)

%c_eval('ePitch?_!bins = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
%c_eval('eFlux_fov?_!bins = ePDist?.tlim(tint_zoom).flux.pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)

%%
npanels = 4;
[h,h2] = initialize_combined_plot(npanels,1,1,0.6,'vertical'); % horizontal

hca = irf_panel('phase space density');
irf_plot(hca,{ePitch1_1bins.elim(500),ePitch1_2bins.elim(500),ePitch1_3bins.elim(500),ePitch1_4bins.elim(500)},'comp')
hca.YLabel.String = {'phase space density',sprintf('(%s)',ePitch1_1bins.units)};


hca = irf_panel('field-of-view flux');
irf_plot(hca,{eFlux_fov1_1bins.elim(500),eFlux_fov1_2bins.elim(500),eFlux_fov1_3bins.elim(500),eFlux_fov1_4bins.elim(500)},'comp')
hca.YLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};

hca = irf_panel('reduced dist');
irf_plot(hca,{eDred1_,eDred2_,eDred3_,eDred4_},'comp')
hca.YLabel.String = {'reduced','phase space density',sprintf('(%s)',eDred1_.units)};

hca = irf_panel('reduced flux');
irf_plot(hca,{eFlux_red1,eFlux_red2,eFlux_red3,eFlux_red4},'comp')
hca.YLabel.String = {'flux of reduced distribution',sprintf('(%s)',eFlux_red1.units)};

irf_zoom(h,'x',tint_zoom)

isub = 1;
hca = h2(isub); isub = isub + 1;
scatter(hca,eFlux_fov1_1bins.elim(500).data,eFlux_red1.data)
if 0
  hold(hca,'on')
  scatter(hca,eFlux_fov1_2bins.elim(500).data,eFlux_red1.data)
  scatter(hca,eFlux_fov1_3bins.elim(500).data,eFlux_red1.data)
  scatter(hca,eFlux_fov1_4bins.elim(500).data,eFlux_red1.data)
  hold(hca,'off')
end

hca.XLabel.String = {'field-of-view flux',sprintf('(%s)',eFlux_fov1_1bins.units)};
hca.YLabel.String = {'flux of reduced distribution',sprintf('(%s)',eFlux_red1.units)};
axis(hca,'equal','square')
hca.XLim = [0 3]*1e6;
hca.YLim = [0 3]*1e6;
hca.Box = 'on';
hca.XGrid = 'on';
hca.YGrid = 'on';

%% Skymaps of FPI, to see if there are significant non-gyrotropies that may explain discrepancies
tint_zoom = irf.tint('2017-07-06T13:54:05.000Z/2017-07-06T13:54:06.000Z');
tint_zoom_edi = irf.tint('2017-07-06T13:54:05.510Z/2017-07-06T13:54:05.630Z'); % 4 ePDist
nrows = 4;
ncols = 4;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

elim = 500;

% Multiple spacecraft
isub = 1;
ic_all = 1:4; % spacecraft id
for ic = ic_all
  c_eval('pdist = ePDist?.tlim(tint_zoom_edi).elim(elim);',ic)
  nt = pdist.length;  
  for it = 1:nt    
    hca = h(isub); isub = isub + 1;
    B = mean(dmpaB1.tlim(pdist(it).time+0.033*0.5*[-1 1]).data,1);
    mms.plot_skymap(hca,pdist(it).dpflux,'log','vectors',{B,'B'});
  end
end
hb = findobj(gcf,'type','colorbar');
delete(hb(2:end))
hb = findobj(gcf,'type','colorbar');
hb.Position(1) = hb.Position(1)+0.03;