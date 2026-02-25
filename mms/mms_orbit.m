%localuser = datastore('local','user');
localuser = 'cecilia';
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   
units = irf_units;

%tint = irf.tint('2015-10-01T00:00:00.00Z/2022-11-05T00:00:00.00Z');
%tic; c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1); toc
%plot(gseR1.data(:,1),gseR1.data(:,2),'-')
%disp('Done.')

hca = subplot(1,1,1);
for year = 2015:2022
  for month = 1:12
    tint = irf.tint(sprintf('%g-%02.0f-01T00:00:00.00Z/%g-%02.0f-05T00:00:00.00Z',year,month,year,month));
    tint
    ts = mms.get_data('R_gse',tint,1);
    if isempty(ts), continue; end
    plot(hca,ts.data(:,1)*1e3/units.RE,ts.data(:,2)*1e3/units.RE)
    hca.Title.String = sprintf('%g-%02.0f',year,month);
    hca.YLim = [-25 25];
    hca.XLim = [-25 25];
    drawnow
    pause(0.1)
  end
end