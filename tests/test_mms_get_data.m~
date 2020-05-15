tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

vars = mms.get_data('vars',tint,1);
nvars = numel(vars);
iTest = find(contains(vars,'Omnifluxion_epd_eis_brst_l2'));

for ivar = iTest
  %ivar
  varstr = vars{ivar}
  eval([varstr ' = mms.get_data(''' varstr ''',tint,1)'])  
end


