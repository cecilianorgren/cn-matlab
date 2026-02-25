

c_eval('tmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l1b/dis-dist/2015/10/16/mms?_fpi_brst_l1b_dis-dist_20151016103254_v1.1.0.cdf'']);',ic);
%c_eval('tmpDataObj? = dataobj(''data/mms?_fpi_brst_l1b_des-dist_20151202011414_v1.1.0.cdf'');',ic);
c_eval('pdist? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_dist''));',ic);
c_eval('energy0? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_energy0'');',ic);
c_eval('energy1? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_energy1'');',ic);
c_eval('phi? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_phi''));',ic);
c_eval('theta? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_theta'');',ic);
c_eval('stepTable? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_stepTable_parity''));',ic);

c_eval('imoments? = mms.psd_moments(pdist?,phi?,theta?,stepTable?,energy0?,energy1?,P?brst,''ion'');',ic);
%% Plot density
h = irf_plot(7);
isub = 1;
for ii = 1:4;
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{ni?brst,imoments?.n_psd,ne?brst,emoments?.n_psd},''comp'')',ii);
  irf_legend(hca,{'i-fpi','i-irfu','e-fpi','e-irfu'},[0.99 0.99])
  ylabel(hca,irf_ssub('n (cm^{-3}) MMS?',ii))
end

hca = h(isub); isub = isub + 1;
irf_plot(hca,{emoments1.n_psd,emoments2.n_psd,emoments3.n_psd,emoments4.n_psd},'comp')
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.99 0.99])
ylabel(hca,'n_e (V) irfu')

hca = h(isub); isub = isub + 1;
irf_plot(hca,'ne?brst','comp')
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.99 0.99])
ylabel(hca,'n_e (V) fpi')

hca = h(isub); isub = isub + 1;
irf_plot(hca,'P?brst','comp')
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.99 0.99])
ylabel(hca,'scPot (V)')