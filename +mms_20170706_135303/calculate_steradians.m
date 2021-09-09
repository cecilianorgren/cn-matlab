% calculate_steradians

% How many steradians is a given pitch angle range

% A = int(sin(theta)*dtheta*dphi)
% dphi is always 2pi

th1 = 0; 
th2 = 11.25*1;

th1 = th1*(pi/180);
th2 = th2*(pi/180);

theta = linspace(th1,th2,1000);
dtheta = theta(2)-theta(1);

f = sin(theta);
intf = sum(f)*dtheta;

str = intf*2*pi