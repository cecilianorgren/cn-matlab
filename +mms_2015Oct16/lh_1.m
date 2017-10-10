disp('Loading magnetic field...')
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint); dmpaB?brst = irf.ts_vec_xyz(dmpaB?brst.time,dmpaB?brst.data(:,1:3)); dmpaB?brst.units = ''nT'';',1:4);
c_eval('dmpaB?brst.name = ''B? brst DMPA'';')

scmPath = '/Users/Cecilia/Data/MMS/2015Oct16/scm/';
scmFiles = dir([scmPath '*.mat']);
for ii = 1:numel(scmFiles)
  load([scmPath scmFiles(ii).name]);
  c_eval('dmpaB?scm = Bscm; dmpaB?scm.units = ''nT''; dmpaB?scm.name = ''scm B?''; dmpaB?scm.coordinateSystem = ''DMPA'';',ii)
end

disp('Loading electric field...')
c_eval('dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('dslE?3Dbrst=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

