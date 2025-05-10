% Partial moments in the tail
localuser = 'cecilia';
localuser = 'cecilianorgren';
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
mms.db_init('local_file_db',['/Volumes/mms']);
ic = 3;
%load /Users/cecilia/'IRFU Dropbox'/'Cecilia Norgren'/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv
%tints = mms_bbfs_db_2017_2022;
%load('/Users/cecilia/IRFU\ Dropbox/Cecilia\ Norgren/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv')

%fid = fopen("/Users/cecilia/IRFU Dropbox/Cecilia Norgren/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv",'r');
%data = fread(fid);

%%
tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',1:4);

c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,ic);',ic)
c_eval('partNi? = mms.get_data(''partNi_fpi_brst_l2'',tint,?);',ic)
c_eval('partVi? = mms.get_data(''partVi_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('partPi? = mms.get_data(''partPi_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('partTi? = mms.get_data(''partTi_gse_fpi_brst_l2'',tint,?);',ic)



%%
nPanels = 4;
h = irf_plot(nPanels);

hca = irf_panel('B');
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'')',ic)
hca.YLabel.String = 'B (nT)';

if 1 % density
  hca = irf_panel('n');
  c_eval('partNi = partNi?;',ic)
  
  
end