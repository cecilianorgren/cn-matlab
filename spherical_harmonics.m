% Define constants.
degree = 2; % l
order = 1; % m

% Create the grid
delta = pi/40;
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);

% Calculate the harmonic
Ymn = legendre(degree,cos(theta(:,1)));

Ymn = Ymn(abs(order)+1,:)';
yy = Ymn;
for kk = 2: size(theta,1)
    yy = [yy Ymn];
end;
yy = yy.*cos(order*phi);

order = max(max(abs(yy)));
rho = ones(size(yy));%60 + 2*yy/order;

for k=0:100;
% Apply spherical coordinate equations
r = rho.*sin(theta);
x = r.*cos(phi+0.1*k);    % spherical coordinate equations
y = r.*sin(phi+0.1*k);
z = rho.*cos(theta);

% Plot the surface
clf
surf(x,y,z,yy)
light
lighting phong
axis tight equal off
view(40,30)
colorbar 

%camzoom(1.5)
pause(0.02)
end