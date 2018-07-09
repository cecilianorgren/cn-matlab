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
ic = 1:4;
Tint = tint;
tmp = ['tmpDataObj? = dataobj(mms.get_filepath(''mms?_edi_brst_l2_amb-pm2'',Tint(1)+20));'];
c_eval(tmp,ic);

% to browse all data provided
c_eval('var = get_variable(tmpDataObj?,''mms?_edi_flux!_0_brst_l2'');',1,1)

edi_nodes = 1:4;
c_eval('flux?_0_! = mms.variable2ts(get_variable(tmpDataObj?,''mms?_edi_flux!_0_brst_l2''));',ic,edi_nodes);
%c_eval('flux0b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_0_brst_l2''));',ic);
c_eval('flux?_180_! = mms.variable2ts(get_variable(tmpDataObj?,''mms?_edi_flux!_180_brst_l2''));',ic,edi_nodes);
%c_eval('flux180b = mms.variable2ts(get_variable(tmpDataObj,''mms?_edi_flux3_180_brst_l2''));',ic);

%flux0180 = irf.ts_scalar(flux180a.time,[(flux0a.data+flux0b.data)/2 (flux180a.data+flux180b.data)/2]);

c_eval('flux0_mms? = [];',ic)
c_eval('flux0_mms? = [flux0_mms? flux?_0_!.data];',ic,edi_nodes)
c_eval('flux0_mms? = irf.ts_scalar(flux?_0_1.time,flux0_mms?);',ic)
c_eval('flux0_mms?.units = var.UNITS; flux0_mms?.siConversion = var.SI_CONVERSION;',ic)

c_eval('flux180_mms? = [];',ic)
c_eval('flux180_mms? = [flux180_mms? flux?_180_!.data];',ic,edi_nodes)
c_eval('flux180_mms? = irf.ts_scalar(flux?_180_1.time,flux180_mms?);',ic)
c_eval('flux180_mms?.units = var.UNITS; flux180_mms?.siConversion = var.SI_CONVERSION;',ic)
%all_flux180 = irf.ts_scalar(flux180_1.time,all_flux180);

c_eval('energy1_mms? = get_variable(tmpDataObj?,''mms?_edi_energy_gdu1_brst_l2'');',ic);
c_eval('energy2_mms? = get_variable(tmpDataObj?,''mms?_edi_energy_gdu2_brst_l2'');',ic);

EDIen = single(energy1_mms1.data(1));