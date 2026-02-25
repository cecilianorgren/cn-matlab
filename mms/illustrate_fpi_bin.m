%clc; clear; close all;

%% --- Bin limits (degrees) ---
theta_min = 22.5;
theta_max = 33.75;

theta_min = 90-11.25*4;
theta_max = theta_min + 11.25;

phi_min = 0;
phi_max = 11.25;

r_min = 1;     % inner radius (energy min)
r_max = 1.3;     % outer radius (energy max)

% Convert to radians
theta_min = deg2rad(theta_min);
theta_max = deg2rad(theta_max);
phi_min   = deg2rad(phi_min);
phi_max   = deg2rad(phi_max);

%% Resolution (smoothness of surfaces)
n_theta = 10;
n_phi   = 10;
n_r     = 2;

theta = linspace(theta_min, theta_max, n_theta);
phi   = linspace(phi_min, phi_max, n_phi);
r     = linspace(r_min, r_max, n_r);

[TH, PH] = meshgrid(theta, phi);

figure; hold on;

%% ---- Outer spherical surface (r_max) ----
R = r_max;
X = R .* sin(TH) .* cos(PH);
Y = R .* sin(TH) .* sin(PH);
Z = R .* cos(TH);

surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

%% ---- Inner spherical surface (r_min) ----
R = r_min;
X = R .* sin(TH) .* cos(PH);
Y = R .* sin(TH) .* sin(PH);
Z = R .* cos(TH);

surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

%% ---- Theta boundary surfaces ----
[RR, PHI] = meshgrid(r, phi);

% theta_min surface
THETA = theta_min * ones(size(RR));
X = RR .* sin(THETA) .* cos(PHI);
Y = RR .* sin(THETA) .* sin(PHI);
Z = RR .* cos(THETA);
surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

% theta_max surface
THETA = theta_max * ones(size(RR));
X = RR .* sin(THETA) .* cos(PHI);
Y = RR .* sin(THETA) .* sin(PHI);
Z = RR .* cos(THETA);
surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

%% ---- Phi boundary surfaces ----
[RR, THETA] = meshgrid(r, theta);

% phi_min surface
PHI = phi_min * ones(size(RR));
X = RR .* sin(THETA) .* cos(PHI);
Y = RR .* sin(THETA) .* sin(PHI);
Z = RR .* cos(THETA);
surf(X,Y,Z,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

% phi_max surface
PHI = phi_max * ones(size(RR));
X = RR .* sin(THETA) .* cos(PHI);
Y = RR .* sin(THETA) .* sin(PHI);
Z = RR .* cos(THETA);
surf(X,Y,Z,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

%% Add random particles
nP = 100;
az = theta_min + (theta_max-theta_min)*rand(nP,1);
pol = phi_min + (phi_max-phi_min)*rand(nP,1);
v = r_min + (r_max-r_min)*rand(nP,1);

x = v .* sin(az) .* cos(pol);
y = v .* sin(az) .* sin(pol);
z = v .* cos(az);
scatter3(x,y,z,'k','filled')
%% ---- Plot formatting ----
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3)
camlight; lighting gouraud
grid on
title('Single Spherical Bin')
hca=gca; 
hca.Visible = 'off';
