% Partial moments in the tail
localuser = 'cecilia';
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
ic = 3;
%load /Users/cecilia/'IRFU Dropbox'/'Cecilia Norgren'/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv
%tints = mms_bbfs_db_2017_2022;
%load('/Users/cecilia/IRFU\ Dropbox/Cecilia\ Norgren/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv')

%fid = fopen("/Users/cecilia/IRFU Dropbox/Cecilia Norgren/Data/Databases/Richard_BBFs/mms_bbfs_db_2017-2022.csv",'r');
%data = fread(fid);
%%
c_eval('partNi? = mms.get_data(''partNi_fpi_fast_l2'',tint,?);',ic)
c_eval('partVi? = mms.get_data(''partVi_fpi_fast_l2'',tint,?);',ic)
c_eval('partPi? = mms.get_data(''partPi_fpi_fast_l2'',tint,?);',ic)
c_eval('partTi? = mms.get_data(''partTi_fpi_fast_l2'',tint,?);',ic)



%%