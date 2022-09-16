tint = [EpochTT('2017-07-25T22:04:44.822634526Z') EpochTT('2017-07-25T22:11:31.820198151Z')];

ic = 2;
c_eval('eis_omni? = mms.get_data(''Omnifluxproton_epd_eis_brst_l2'',tint,?);',2)
c_eval('eis_pa? = mms.get_data(''Pitchanglefluxproton_epd_eis_brst_l2'',tint,?);',2)
c_eval('feeps_ion_omni? = mms.get_data(''Omnifluxion_epd_feeps_brst_l2'',tint,?);',1:4)
c_eval('feeps_ion_pa? = mms.get_data(''Pitchanglefluxion_epd_feeps_brst_l2'',tint,?);',1:4)
%% 4 sc average (3 in case of eis since mms1 seems to be missing data)
% Feeps
omnidata1234 = (feeps_ion_omni1.data + ...
                feeps_ion_omni2.data(1:end-1,:) + ...
                feeps_ion_omni3.data + ...
                feeps_ion_omni4.data)/4;
feeps_ion_omni1234 =  PDist(feeps_ion_omni1.time,omnidata1234,'omni',feeps_ion_omni1.depend{1}); % energies keV -> eV            
feeps_ion_omni1234.siConversion = feeps_ion_omni1.siConversion;
feeps_ion_omni1234.units = feeps_ion_omni1.units;
feeps_ion_omni1234.species = feeps_ion_omni1.species;

nPointTotal = feeps_ion_pa1.ancillary.Coverage.N + ...
  feeps_ion_pa2.ancillary.Coverage.N(1:end-1,:,:) +...
  feeps_ion_pa3.ancillary.Coverage.N + ...
  feeps_ion_pa4.ancillary.Coverage.N;

padata1234 = (feeps_ion_pa1.data.*feeps_ion_pa1.ancillary.Coverage.N + ...
                feeps_ion_pa2.data(1:end-1,:,:).*feeps_ion_pa2.ancillary.Coverage.N(1:end-1,:,:) + ...
                feeps_ion_pa3.data.*feeps_ion_pa3.ancillary.Coverage.N + ...
                feeps_ion_pa4.data.*feeps_ion_pa4.ancillary.Coverage.N)./nPointTotal;
feeps_ion_pa1234 =  PDist(feeps_ion_pa1.time,padata1234,'pitchangle',feeps_ion_pa1.depend{1},feeps_ion_pa1.depend{2}); % energies keV -> eV            
feeps_ion_pa1234.siConversion = feeps_ion_pa1.siConversion;
feeps_ion_pa1234.units = feeps_ion_pa1.units;
feeps_ion_pa1234.species = feeps_ion_pa1.species;

% EIS
omnidata1234 = (...
                eis_omni2.data(1:end-1,:) + ...
                eis_omni3.data + ...
                eis_omni4.data)/4;
eis_omni1234 =  PDist(eis_omni3.time,omnidata1234,'omni',eis_omni3.depend{1}); % energies keV -> eV            
eis_omni1234.siConversion = eis_omni3.siConversion;
eis_omni1234.units = eis_omni3.units;
eis_omni1234.species = eis_omni3.species;

nPointTotal = ...
  eis_pa2.ancillary.Coverage.N(1:end-1,:,:) +...
  eis_pa3.ancillary.Coverage.N + ...
  eis_pa4.ancillary.Coverage.N;

padata1234 = (eis_pa2.data(1:end-1,:,:).*eis_pa2.ancillary.Coverage.N(1:end-1,:,:) + ...
                eis_pa3.data.*eis_pa3.ancillary.Coverage.N + ...
                eis_pa4.data.*eis_pa4.ancillary.Coverage.N)./nPointTotal;
eis_pa1234 =  PDist(eis_pa3.time,padata1234,'pitchangle',eis_pa3.depend{1},eis_pa3.depend{2}); % energies keV -> eV            
eis_pa1234.siConversion = eis_pa3.siConversion;
eis_pa1234.units = eis_pa3.units;
eis_pa1234.species = eis_pa3.species;
%% Plot
h = irf_plot(8);
elim_eis = [32000 Inf];
elim_eis = [0 Inf];
if 1 % FEEPS Pitch all sc
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,eis_pa1234.elim(elim_eis).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 1 % i DEF feeps omni all sc
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,eis_omni1234.elim(elim_eis).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
end
for ic = 2:4
  if 1 % EIS Pitch 0 allE
    isub = isub + 1;
    hca = irf_panel(['eis Pitch all E mms' num2str(ic)]);  
    %elim = [0 32000];
    c_eval('irf_spectrogram(hca,eis_pa?.elim(elim_eis).specrec(''pitchangle''));',ic)
    hca.YLim = [0 180];
    hca.YLabel.String = {'\theta','(deg)'};  
    hca.YLabel.Interpreter = 'tex';
    %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  end
  if 1 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel(['i DEF eis omni mms' num2str(ic)]);  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.elim(elim_eis).specrec(''energy''),''log'');',ic)
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
end
elim_feeps = [8e04 Inf];
if 0 % FEEPS Pitch all sc
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 0 % i DEF feeps omni all sc
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
for ic = []%1:4
  if 1 % FEEPS Pitch 0 allE
  %  feeps_pa2.elim([9.7+04 Inf])
    isub = isub + 1;
    hca = irf_panel(['feeps pitch mms ' num2str(ic)]);  
    %elim = [0 32000];
    c_eval('irf_spectrogram(hca,feeps_ion_pa?.elim(elim_feeps).specrec(''pitchangle''));',ic)
    hca.YLim = [0 180];
    hca.YLabel.String = {'\theta','(deg)'};  
    hca.YLabel.Interpreter = 'tex';
    %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  end
  if 1 % i DEF feeps omni
    isub = isub + 1;
    hca = irf_panel(['feeps omni mms ' num2str(ic)]);  
    c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.elim(elim_feeps).specrec(''energy''),''log'');',ic)
    set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
    colormap(hca,cmap) 
    hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
  end
end

colormap(obs_colors('candy4'))
%hlinks = linkprop(h(3:end),{'CLim'});
hlinks = linkprop(h(1:6),{'CLim'});