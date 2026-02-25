units = irf_units;
B0 = 1;
J0 = 1;
n0 = 1;
v0 = 1;
L = 1;

B = @(z) B0*tanh(z/L); %BofL = @(z,L) B0*tanh(z/L);
n = @(z) n0*cosh(z/L).^-2;
J = @(z) J0*(1 - tanh(z/L).^2);
v = @(z,u) n0*u*cosh(z/L).^-2;
 
z = linspace(-3*L,3*L,100);

%% B n
hca = subplot(1,1,1);
plot(hca,z,B(z),z,B(z).^2,z,n(z),z,J(z),'--');
hca.FontSize = 16;
hca.XLabel.String = 'y/L';
hca.YLabel.String = 'Normalized amplitudes';
hleg = legend('B','|B|','n','J');
hleg.Box = 'off';
%% v
colors = mms_colors('matlab');
Ti = 1;
Te = 1;
ui = 1;
ue = -ui*Te/Ti;

hca = subplot(1,1,1);
hca.ColorOrder = colors(1:2,:);
hca.LineStyleOrder = {'-','--'};
lines = plot(hca,z,v(z,ui),z,v(z,ue),z,v(z,0.5*ue));
lines(3).Color = colors(2,:);
lines(3).LineStyle = '--';
hca.FontSize = 16;
hca.XLabel.String = 'y/L';
hca.YLabel.String = 'Normalized velocities';
hleg = legend('v_i','v_e: T_e = T_i','v_e: T_e = 0.5 T_i');
hleg.Box = 'off';

%% B n v
hca = subplot(2,1,1);
plot(hca,z,B(z),z,B(z).^2,z,n(z));
hca.FontSize = 16;
hca.XLabel.String = 'y/L';
hca.YLabel.String = 'Normalized amplitudes';
hleg = legend('B','|B|','n');
hleg.Box = 'off';

hca = subplot(2,1,2);
hca.ColorOrder = colors(1:2,:);
hca.LineStyleOrder = {'-','--'};
lines = plot(hca,z,v(z,ui),z,v(z,ue),z,v(z,0.5*ue));
lines(3).Color = colors(2,:);
lines(3).LineStyle = '--';
hca.FontSize = 16;
hca.XLabel.String = 'y/L';
hca.YLabel.String = 'Normalized velocities';
hleg = legend('v_i','v_e: T_e = T_i','v_e: T_e = 0.5 T_i');
hleg.Box = 'off';