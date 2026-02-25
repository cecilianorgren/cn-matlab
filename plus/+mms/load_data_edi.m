% tmp = ['tmpDataObj = dataobj(mms.get_filepath(''mms?_edi_brst_l2_amb-pm2'',Tint));'];
% c_eval(tmp,ic);
% 
% c_eval('flux0a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_0_brst_l2''));',ic);
% c_eval('flux0b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_0_brst_l2''));',ic);
% c_eval('flux180a = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux2_180_brst_l2''));',ic);
% c_eval('flux180b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_180_brst_l2''));',ic);
% 
% flux0180 = irf.ts_scalar(flux180a.time,[(flux0a.data+flux0b.data)/2 (flux180a.data+flux180b.data)/2]);
% 
% c_eval('energy1 = get_variable(tmpDataObj,''mms?_edi_energy_gdu1_brst_l2'');',ic);
% c_eval('energy2 = get_variable(tmpDataObj,''mms?_edi_energy_gdu2_brst_l2'');',ic);
% 
% EDIen = single(energy1.data(1));

%%
%ic = 1:4;
units = irf_units;
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
pa_edi_edges = 0:11.25:45;
pa_edi_centers = pa_edi_edges(1:4) + 11.25/2;
v_edi_edges = [v_edi_minus v_edi_plus];
v_edi_par_edges = tocolumn(v_edi_edges)*torow(cosd(pa_edi_edges));
v_edi_par = v_edi*cosd(pa_edi_centers);
v_edi_par_edges_center = v_edi*cosd(pa_edi_edges);


Tint = tint;

tmp = ['tmpDataObj? = dataobj(mms.get_filepath(''mms?_edi_brst_l2_amb-pm2'',Tint(1)+20));'];
%tmp = ['tmpDataObj? = dataobj(mms.get_filepath(''mms?_edi_brst_l2_amb'',Tint(1)+20));'];

c_eval(tmp,ic);

% to browse all data provided
c_eval('var = get_variable(tmpDataObj?,''mms?_edi_flux!_0_brst_l2'');',1,1)

edi_nodes = 1:4;
c_eval('flux?_0_! = mms.variable2ts(get_variable(tmpDataObj?,''mms?_edi_flux!_0_brst_l2''));',ic,edi_nodes);
%c_eval('flux0b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_0_brst_l2''));',ic);
c_eval('flux?_180_! = mms.variable2ts(get_variable(tmpDataObj?,''mms?_edi_flux!_180_brst_l2''));',ic,edi_nodes);
%c_eval('flux180b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_180_brst_l2''));',ic);

pa_edi = 0:11.24:45; % scaling the same as for 0-45 as for 180-135
pa_scaling = pi*(cosd(pa_edi(1:4)).^2-cosd(pa_edi(2:5)).^2);  

c_eval('flux_0_45_mms? = irf.ts_scalar(flux?_0_1.time,[flux?_0_1.data flux?_0_2.data flux?_0_3.data flux?_0_4.data]); flux_0_45_mms?.units =''cm^-2 s^-1 sr^-1'';',ic)
c_eval('flux_0_45_scaled_mms? = irf.ts_scalar(flux_0_45_mms?.time,flux_0_45_mms?.data.*repmat(pa_scaling,flux_0_45_mms?.length,1)); flux_0_45_scaled_mms?.units = ''cm^-2 s^-1'';',ic)
c_eval('flux_par_int_mms? = irf.ts_scalar(flux_0_45_scaled_mms?.time,sum(flux_0_45_scaled_mms?.data,2)); flux_par_int_mms?.units = ''cm^-2 s^-1'';',ic)

c_eval('flux_135_180_mms? = irf.ts_scalar(flux?_180_1.time,[flux?_180_1.data flux?_180_2.data flux?_180_3.data flux?_180_4.data]); flux_135_180_mms?.units =''cm^-2 s^-1 sr^-1'';',ic)
c_eval('flux_135_180_scaled_mms? = irf.ts_scalar(flux_135_180_mms?.time,flux_135_180_mms?.data.*repmat(pa_scaling,flux_135_180_mms?.length,1)); flux_135_180_scaled_mms?.units = ''cm^-2 s^-1'';',ic)
c_eval('flux_apar_int_mms? = irf.ts_scalar(flux_135_180_scaled_mms?.time,sum(flux_135_180_scaled_mms?.data,2)); flux_apar_int_mms?.units = ''cm^-2 s^-1'';',ic)


%flux0180 = irf.ts_scalar(flux180a.time,[(flux0a.data+flux0b.data)/2 (flux180a.data+flux180b.data)/2]);

c_eval('flux0_mms? = [];',ic)
c_eval('flux0_mms? = [flux0_mms? flux?_0_!.data];',ic,edi_nodes)
c_eval('flux0_mms? = irf.ts_scalar(flux?_0_1.time,flux0_mms?);',ic)
c_eval('flux0_mms?.units = var.UNITS; flux0_mms?.siConversion = var.SI_CONVERSION;',ic)
%c_eval('flux0_mms?.ancillary.nodes = edi_nodes;',ic)

c_eval('flux180_mms? = [];',ic)
c_eval('flux180_mms? = [flux180_mms? flux?_180_!.data];',ic,edi_nodes)
c_eval('flux180_mms? = irf.ts_scalar(flux?_180_1.time,flux180_mms?);',ic)
c_eval('flux180_mms?.units = var.UNITS; flux180_mms?.siConversion = var.SI_CONVERSION;',ic)
%c_eval('flux180_mms?.ancillary.nodes = edi_nodes;',ic)
%all_flux180 = irf.ts_scalar(flux180_1.time,all_flux180);

c_eval('energy1_mms? = get_variable(tmpDataObj?,''mms?_edi_energy_gdu1_brst_l2'');',ic);
c_eval('energy2_mms? = get_variable(tmpDataObj?,''mms?_edi_energy_gdu2_brst_l2'');',ic);

EDIen = single(energy1_mms1.data(1));