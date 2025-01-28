tint = irf.tint('2015-10-16T10:30:00.00Z/2015-10-16T10:30:40.00Z');
tint = irf.tint('2015-10-16T10:30:27.00Z/2015-10-16T10:30:27.50Z');

% load particle distributions
tic; c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ic); toc
tic; c_eval('disti = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint);',ic); toc

%%
% particle distribution angles
dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;

% particle distribution energies
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);



% define directions