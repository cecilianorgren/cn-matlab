%% 
units = irf_units;
pos = load('/Users/cecilianorgren/Data/Databases/DB_Lalti/Proba_full_mms1_pos.mat');
gseR = pos.dTmpR;
time = EpochTT('2017-07-16T12:00:00');

r = gseR.tlim(time+60*60*1*[-4 0]);


plot(r.x.data*1e3/units.RE,r.y.data*1e3/units.RE)
axis equal

% '2017-07-16T02:09:14.661335000Z'
% '2017-07-16T03:06:10.186133000Z'
% '2017-07-16T06:25:31.272540000Z'
% '2017-07-16T06:27:10.273254000Z'
% '2017-07-16T06:33:23.775943000Z'
% '2017-07-16T06:35:11.776722000Z'
% '2017-07-16T06:45:01.280985000Z'