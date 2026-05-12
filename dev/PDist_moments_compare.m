mms.db_init('local_file_db','/Volumes/mms'); % set up database
units = irf_units;
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); % Torbert event 
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)

c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);

%% Remove noise
nMean = [5,3,3,3]; nThresh = 3;

c_eval('PD_counts = iPDist?_counts;',ic)
c_eval('PD_orig = iPDist?;',ic)

nCounts = 5;
nMovMean = [7 7];
matEmask = PD_counts.find_low_counts('counts',nCounts,'nMovMean',nMovMean,'output','mat');

PD_notmasked = PD_orig;

PD_masked = PD_orig.mask('energy','mat',matEmask);

h = irf_plot({PD_orig.deflux.omni.specrec,PD_masked.deflux.omni.specrec}); c_eval('h(?).YScale = ''log'';',1:numel(h))
hlinks = linkprop(h,{'CLim','YLim'});



%%
% Calculate moments
n_orig = PD_orig.n;
n = PD_masked.n;

v_orig = PD_orig.vel;
v = PD_masked.vel;

p_orig = PD_orig.p;
p = PD_masked.p;

t_orig = (p_orig*1e-9/(n_orig*1e6))/units.eV;
t = (p*1e-9/(n*1e6))/units.eV;
%%
h = irf_plot(6);

hca = irf_panel('DEF omni orig');
irf_spectrogram(hca,PD_orig.deflux.omni.specrec)
irf_legend(hca,{'orig VDF'},[0.98 0.05],'color','k')

hca = irf_panel('DEF omni masked');
irf_spectrogram(hca,PD_masked.deflux.omni.specrec)
irf_legend(hca,{'clean VDF'},[0.98 0.05],'color','k')

hca = irf_panel('n');
irf_plot(hca,{n_orig,n,ni3,ne3},'comp')
hca.YLabel.String = 'n (cc)';
irf_legend(hca,{'orig VDF','clean VDF','n_i FPI','n_e FPI'},[0.98 0.98])

hca = irf_panel('v');
irf_plot(hca,{v_orig.x,v.x,gseVi3.x},'comp')
hca.YLabel.String = 'v_x (km/s)';
irf_legend(hca,{'orig VDF','clean VDF','FPI'},[0.98 0.98])

hca = irf_panel('P');
irf_plot(hca,{p_orig.trace/3,p.trace/3,gsePi3.trace/3},'comp')
hca.YLabel.String = 'P (nPa)';
irf_legend(hca,{'orig VDF','clean VDF','FPI'},[0.98 0.98])

hca = irf_panel('T');
irf_plot(hca,{t_orig.trace/3,t.trace/3,gseTi3.trace/3},'comp')
hca.YLabel.String = 'T (eV)';
irf_legend(hca,{'orig VDF','clean VDF','FPI'},[0.98 0.98])

irf_plot_axis_align
irf_zoom(h,'x',tint)
h(end).XTickLabelRotation = 0;