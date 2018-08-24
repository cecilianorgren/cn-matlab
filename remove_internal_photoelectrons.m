% bgdist model is given with start_del_phi as depend_0. That is, bgdist
% changes with the phase of the spacecraft (because the spacecraft is not 
% spherically symmetric).
% 'Nominally, spin-phase counts range over [0000, 5759] cnts.  Each count 
% represents 1/16 deg of observatory spin-phase.  Zero spin-phase count 
% indicates obs +X-axis aligned with Sun.'
% 5759/16 = 359.9375
% So we must loop through a ePDist, find the phase, and assign the correct
% bgdist based on that: start_del_phi -> time

%% Load FPI data to compare to
tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
mms_id = 1;
c_eval('pathPDist = ''/Users/cecilia/Data/fpi/mms?_fpi_brst_l2_des-dist_20170706135303_v3.2.0.cdf'';',mms_id)
ePDist = mms.make_pdist(pathPDist);
ePDist = ePDist.tlim(tint);
%c_eval('ePDist = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

c_eval('start_del_phi = get_variable(dataobj(pathPDist),''mms?_des_startdelphi_count_brst'');',mms_id)
ts_start_del_phi = mms.variable2ts(start_del_phi);
ts_start_del_phi = ts_start_del_phi.tlim(tint);

%% Load background models
datapath = '/Users/cecilia/Data/fpi_bgdist/';
files = dir([datapath '*.cdf']);
imodel = 1;
file = [files(imodel).folder filesep files(imodel).name];
tmpDataObj = dataobj(file);
% fileInfo = regexp(tmpDataObj.GlobalAttributes.Logical_source{1}, ...
%   'mms_fpi_(?<tmMode>(fast|brst))_(?<datalevel>\w*)_(?<detector>\w*)', 'names');

fileInfo = regexp(files(imodel).name, ...
  'mms_fpi_(?<tmMode>(fast|brst))_(?<datalevel>\w*)_(?<detector>\w*)-bgdist_(?<version>\w*.\w*.\w*)_(?<model>\w*-\w*).cdf', 'names');

bgdist_p0 = get_variable(tmpDataObj,'mms_des_bgdist_p0_brst');
bgdist_p1 = get_variable(tmpDataObj,'mms_des_bgdist_p1_brst');

%% Assign right phase bg dist
bgdist = bgdist_p0;
%BGDist_p0 = BGDist;
BGDist = ePDist;
vec_del_phi = double(bgdist.DEPEND_0.data);

%mat_del_phi = repmat(vec_del_phi);
%ind_del_phi = find(abs(vec_del_phi-tmp_del_phi) == min(abs(vec_del_phi-tmp_del_phi)));
fprintf('it = %4.0f/%4.0f\n',0,ePDist.length) % display progress
for it = 1:ePDist.length
  if mod(it,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],it,ePDist.length); end % display progress
  tmp_del_phi = double(ts_start_del_phi.data(it));
  ind_del_phi = find(abs(vec_del_phi-tmp_del_phi) == min(abs(vec_del_phi-tmp_del_phi)));
  BGDist.data(it,:,:,:) = bgdist.data(ind_del_phi(1),:,:,:);
end

pdist_diff = ePDist;
pdist_diff.data = ePDist.data - BGDist.data;
pdist_diff.data(pdist_diff.data<0) = 0;

%% Make reduced distributions
eint = [000 40000];
vint = [-Inf Inf];

scpot = scPot1.resample(ePDist);
lowerelim = scpot/scpot*0;
eLine = dmpaB1.resample(ePDist).norm;
tic; 
ef1D = ePDist.reduce('1D',dmpaB1.resample(ePDist).norm,'vint',vint,'scpot',scPot1.resample(ePDist));
ef1D_nobg = pdist_diff.reduce('1D',dmpaB1.resample(pdist_diff).norm,'vint',vint,'scpot',scPot1.resample(pdist_diff));
ef1D_bg = BGDist.reduce('1D',dmpaB1.resample(BGDist).norm,'vint',vint,'scpot',scPot1.resample(BGDist));
toc % reduced distribution along B  

%% plot
h = irf_plot(2);

if 0 % ePDist omni
  hca = irf_panel('epdist omni');
  irf_spectrogram(hca,ePDist.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{'original'},[0.02 0.98],'color',[1 1 1])
end
if 0 % BGDist omni
  hca = irf_panel('epdist bg omni');
  irf_spectrogram(hca,BGDist.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{ sprintf('internal photoelectrons, model %s',fileInfo.model);},[0.02 0.98],'color',[1 1 1])
end
if 0 % ePDist-BGDist omni
  hca = irf_panel('epdist omni,r bg subt.');
  irf_spectrogram(hca,pdist_diff.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{'original - internal photoelectrons model'},[0.02 0.98],'color',[1 1 1])
end
if 0 % ePDist omni deflux
  hca = irf_panel('epdist omni');
  irf_spectrogram(hca,ePDist.deflux.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{'original'},[0.02 0.98],'color',[1 1 1])
end
if 0 % BGDist omni deflux
  hca = irf_panel('epdist bg omni');
  irf_spectrogram(hca,BGDist.deflux.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{ sprintf('internal photoelectrons, model %s',fileInfo.model);},[0.02 0.98],'color',[1 1 1])
end
if 0 % ePDist-BGDist omni deflux
  hca = irf_panel('epdist omni,r bg subt.');
  irf_spectrogram(hca,pdist_diff.deflux.omni.specrec)
  hca.YScale = 'log';
  if 1 % add scpot
    hold(hca,'on')
    c_eval('irf_plot(hca,scPot?);',mms_id)
    hold(hca,'off')
  end
  hca.YLabel.String = 'E (eV)';
  irf_legend(hca,{'original - internal photoelectrons model'},[0.02 0.98],'color',[1 1 1])
end
if 1 % e psd vpar
  hca = irf_panel('ePDist reduced');
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));    
  hca.YLabel.String = 'v_e (km/s)'; 
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % e psd vpar
  hca = irf_panel('ePDist reduced bg');
  irf_spectrogram(hca,ef1D_bg.specrec('velocity_1D','10^3 km/s'));    
  hca.YLabel.String = 'v_e (km/s)'; 
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % e psd vpar
  hca = irf_panel('ePDist reduced dist-bg');
  irf_spectrogram(hca,ef1D_nobg.specrec('velocity_1D','10^3 km/s'));    
  hca.YLabel.String = 'v_e (km/s)'; 
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end

for ipanel = 1:numel(h)
  %h(ipanel).CLim = [-36 -24]; % psd
  %h(ipanel).CLim = [3 8]; % deflux
  h(ipanel).CLim = [-8 -2]; % deflux
end

h(1).Title.String = sprintf('BGDist model: %s',fileInfo.model);

%%

bgdist_p0 = get_variable(tmpDataObj,'mms_des_bgdist_p0_brst');
bgdist_p1 = get_variable(tmpDataObj,'mms_des_bgdist_p1_brst');

bulkx_p0 = get_variable(tmpDataObj,'mms_des_bulkx_dbcs_p0_brst');
bulkx_p1 = get_variable(tmpDataObj,'mms_des_bulkx_dbcs_p1_brst');
%%
%tmpDist = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_dist_' fileInfo.tmMode]);
bgdist_p0 = mms_des_bgdist_p0_brst;
bgdist_p1 = mms_des_bgdist_p1_brst;


% Collect User data
ud = [];
ud.GlobalAttributes = tmpDist.GlobalAttributes;
ud.CATDESC          = tmpDist.CATDESC;
if isfield(tmpDist,'DISPLAY_TYPE'), ud.DISPLAY_TYPE     = tmpDist.DISPLAY_TYPE; end
ud.FIELDNAM         = tmpDist.FIELDNAM;
ud.VALIDMIN         = tmpDist.VALIDMIN;
ud.VALIDMAX         = tmpDist.VALIDMAX;
if isfield(tmpDist,'LABLAXIS'), ud.LABLAXIS = tmpDist.LABLAXIS; end
if isfield(tmpDist,'LABL_PTR_1'), ud.LABL_PTR_1 = tmpDist.LABL_PTR_1;
elseif isfield(tmpDist,'LABL_PTR_2'), ud.LABL_PTR_2 = tmpDist.LABL_PTR_2;
elseif isfield(tmpDist,'LABL_PTR_3'), ud.LABL_PTR_3 = tmpDist.LABL_PTR_3;
end