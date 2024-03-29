mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
tint = irf.tint('2015-11-30T00:00:00.00Z/2015-11-30T02:00:00.00Z'); % Eriksson2016 (Elin)
ic = 1;

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_fast_l2'',''mms?_fgm_b_gse_fast_l2'',tint);',ic);
c_eval('gseE?=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_gse_fast_l2'',tint);',ic);

% Plasma distributions
c_eval('iPDist? = mms.get_data(''PDi_fpi_fast_l2'',tint,?);',ic)

% Pressure and temperature
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_fast_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_fast_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_fast_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_fast_l2'',tint,?);',ic)

%c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
%c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
%c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
%c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
c_eval('ne? = mms.get_data(''Ne_fpi_fast_l2'',tint,?);',ic)
c_eval('ni? = mms.get_data(''Ni_fpi_fast_l2'',tint,?);',ic)
c_eval('nhplus? = mms.get_data(''Nhplus_hpca_srvy_l2'',tint,?);',ic)
c_eval('nheplus? = mms.get_data(''Nheplus_hpca_srvy_l2'',tint,?);',ic)
c_eval('noplus? = mms.get_data(''Noplus_hpca_srvy_l2'',tint,?);',ic)


% Velocity
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_fast_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_fast_l2'',tint,?);',ic)

%% Figure
npanels = 2;
h = irf_plot(npanels);
ic = 1;

if 1 % all densities
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1234gb'))
  c_eval('irf_plot(hca,{ne?,ni?,nhplus?,nheplus?,noplus?},''comp'');',ic)
  irf_legend(hca,{'n_e^{FPI}','n_i^{FPI}','n_{H+}^{HPCA}','n_{He+}^{HPCA}','n_{O+}^{HPCA}'},[0.98 0.98])
end

if 1 % hpca densities
  hca = irf_panel('n hpca');
  set(hca,'ColorOrder',mms_colors('1234gb'))
  c_eval('irf_plot(hca,{nhplus?,nheplus?,noplus?},''comp'');',ic)
  irf_legend(hca,{'n_{H+}^{HPCA}','n_{He+}^{HPCA}','n_{O+}^{HPCA}'},[0.98 0.98])
end