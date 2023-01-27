tint = [EpochTT('2017-07-25T22:04:44.822634526Z') EpochTT('2017-07-25T22:11:31.820198151Z')];

ic = 1:4;
c_eval('eis_omni? = mms.get_data(''Omnifluxproton_epd_eis_brst_l2'',tint,?);',ic)
c_eval('eis_pa? = mms.get_data(''Pitchanglefluxproton_epd_eis_brst_l2'',tint,?);',ic)
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


% Different ways to add pitch angle distributions together.
% 'Coverage' specifies which bin had datapoints in them.
% At first I thought it was what bins were were covered by the instrument, 
% i.e. not what what bins have data, but what bins potentially could have data
% (since there is not complete skymap coverage). However, if it was like
% this, a given pitch angle would have coverage for all energies, but this
% is not the case... For example, do
% 
% >> feeps_ion_pa1.ancillary.Coverage.N(1,:,:)
% 
% So how does one do averaging, if one does not know if
% a zero is due to no data in a covered bin, or if there simply was no bin
% at that pitchangle at that time...
%
% feeps_ion_pa1.ancillary.Coverage
%
% ans = 
%
%  struct with fields:
%
%    STR: 'N - Number of data points in each (time,energy,pitchangle) bin.'
%      N: [1323×16×12 double]


% Check what cioverage actually mean. I.e. is there coverage with the data
% being zero?
pa1 = feeps_ion_pa1.data;
cov1 = feeps_ion_pa1.ancillary.Coverage.N;
covered_data = pa1(cov1==1);

% 1
padata1234 =   (feeps_ion_pa1.data.*feeps_ion_pa1.ancillary.Coverage.N + ...
                feeps_ion_pa2.data(1:end-1,:,:).*feeps_ion_pa2.ancillary.Coverage.N(1:end-1,:,:) + ...
                feeps_ion_pa3.data.*feeps_ion_pa3.ancillary.Coverage.N + ...
                feeps_ion_pa4.data.*feeps_ion_pa4.ancillary.Coverage.N)./nPointTotal;
% 2
if 0
padata1234 =   (feeps_ion_pa1.data + ...
                feeps_ion_pa2.data(1:end-1,:,:) + ...
                feeps_ion_pa3.data + ...
                feeps_ion_pa4.data)/4;
end
% 3
if 1
  padata1234 = cat(4,feeps_ion_pa1.data,feeps_ion_pa2.data(1:end-1,:,:),feeps_ion_pa3.data,feeps_ion_pa4.data);  
  padata1234 = nanmean(padata1234,4);
end 
% 4
if 0
  pa1 = feeps_ion_pa1.data;              pa1(feeps_ion_pa1.ancillary.Coverage.N==0) = NaN;  
  pa2 = feeps_ion_pa2.data(1:end-1,:,:); pa2(feeps_ion_pa2.ancillary.Coverage.N(1:end-1,:,:)==0) = NaN;
  pa3 = feeps_ion_pa3.data;              pa3(feeps_ion_pa3.ancillary.Coverage.N==0) = NaN;
  pa4 = feeps_ion_pa4.data;              pa4(feeps_ion_pa4.ancillary.Coverage.N==0) = NaN;
  padata1234 = cat(4,pa1,pa2,pa3,pa4);
  padata1234 = nanmean(padata1234,4);
  %padata1234 = pa1;
end

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


padata1234 = cat(4,eis_pa2.data(1:end-1,:,:),eis_pa3.data,eis_pa4.data);
padata1234 = nanmean(padata1234,4);
eis_pa1234 =  PDist(eis_pa3.time,padata1234,'pitchangle',eis_pa3.depend{1},eis_pa3.depend{2}); % energies keV -> eV            
eis_pa1234.siConversion = eis_pa3.siConversion;
eis_pa1234.units = eis_pa3.units;
eis_pa1234.species = eis_pa3.species;

% Plot
h = irf_plot(8);
elim_eis = [32000 Inf];
elim_eis = [0 Inf];
if 1 % EIS Pitch all sc
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('eis pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,eis_pa1234.elim(elim_eis).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{EIS}','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 1 % i DEF eis omni all sc
  isub = isub + 1;
  hca = irf_panel('eis omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,eis_omni1234.elim(elim_eis).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
end
for ic = 2%2:4
  if 1 % EIS Pitch 0 allE
    isub = isub + 1;
    hca = irf_panel(['eis Pitch all E mms' num2str(ic)]);  
    %elim = [0 32000];
    c_eval('irf_spectrogram(hca,eis_pa?.elim(elim_eis).specrec(''pitchangle''));',ic)
    hca.YLim = [0 180];
    hca.YLabel.String = {'\theta','(deg)'};  
    hca.YLabel.Interpreter = 'tex';
    %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
    irf_legend(hca,sprintf('MMS%g',ic),[0.02 0.98])
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
  irf_legend(hca,sprintf('MMS%g',ic),[0.02 0.98])
end
end
elim_feeps = [8e04 Inf];
if 1 % FEEPS Pitch all sc
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
if 1 % i DEF feeps omni all sc
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
for ic = 2
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
    irf_legend(hca,sprintf('MMS%g',ic),[0.02 0.98])
  end
  if 1 % i DEF feeps omni
    isub = isub + 1;
    hca = irf_panel(['feeps omni mms ' num2str(ic)]);  
    c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.elim(elim_feeps).specrec(''energy''),''log'');',ic)
    set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
    colormap(hca,cmap) 
    hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
    irf_legend(hca,sprintf('MMS%g',ic),[0.02 0.98])
  end
end

colormap(obs_colors('candy4'))
%hlinks = linkprop(h(3:end),{'CLim'});
hlinks = linkprop(h(1:4),{'CLim'});
hlinks = linkprop(h(5:8),{'CLim'});
linkprop(h,{'XLim'});
%linkprop(h(2:2:end),{'YLim'}); h(2).YLim = [eis_omni2.depend{1}(1,end) feeps_ion_omni2.depend{1}(1,end)];
%h(1).CLim = [0 3];