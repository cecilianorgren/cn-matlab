loadPath = '/Users/Cecilia/Data/MMS/2015Oct16/mat-data/';
listMats = dir(loadPath);
%%
nMats = numel(listMats);
for ii=1:nMats
  if ~listMats(ii).isdir && ~strcmp(listMats(ii).name(1),'.')
    load([loadPath listMats(ii).name])
  end
end
%%
scmPath = '/Users/Cecilia/Data/MMS/2015Oct16/scm/';
scmFiles = dir([scmPath '*.mat']);
for ii = 1:numel(scmFiles)
  load([scmPath scmFiles(ii).name]);
  c_eval('dmpaB?scm = Bscm; dmpaB?scm.units = ''nT''; dmpaB?scm.name = ''scm B?''; dmpaB?scm.coordinateSystem = ''DMPA'';',ii)
end
%%
loadPath = '/Users/Cecilia/Data/MMS/2015Oct16/recalculated_moments/';

c_eval('load([loadPath ''psd_emomentsMMS?.mat'']);')
c_eval('Ve_psd? = irf.ts_vec_xyz(psd_emomentsMMS?.Ve_psd?.time,psd_emomentsMMS?.Ve_psd?.data); Ve_psd?.name = ''recalc Ve'';')
c_eval('Pe_psd? = psd_emomentsMMS?.Pe_psd?;')
c_eval('Te_psd? = psd_emomentsMMS?.Te_psd?;')
c_eval('ne_psd? = psd_emomentsMMS?.ne_psd?;')
