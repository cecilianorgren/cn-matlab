ic = 1;
id = 0;

id = id+1;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
c_eval('tic; disDist?_! = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint); toc',ic,id);
c_eval('disDist?_!, disDist?_!.userData.GlobalAttributes.Logical_file_id',ic,id)

id = id+1;
tint = irf.tint('2015-10-16T10:33:00.00Z/2015-10-16T10:34:20.00Z');
c_eval('tic; disDist?_! = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint); toc',ic,id);
c_eval('disDist?_!, disDist?_!.userData.GlobalAttributes.Logical_file_id',ic,id)

id = id+1;

tint = irf.tint('2015-10-16T10:28:00.00Z/2015-10-16T10:40:00.00Z');
c_eval('tic; disDist?_! = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint); toc',ic,id);
c_eval('nDist = numel(disDist?_!);',ic,id)
for ii = 1:nDist
  c_eval('disDist?_!{ii}, disDist?_!{ii}.userData.GlobalAttributes.Logical_file_id',ic,id)
end

%% These should be the same, but are not
disDist1_3{3}.time
disDist1_2.time

disDist1_3
disDist1_3{3}
disDist1_2

