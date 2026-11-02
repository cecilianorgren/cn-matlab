%clc; clear; close all;

%% --- Bin limits (degrees) ---
theta_min = 22.5;
theta_max = 33.75;

theta_min = 90;
theta_max = 102.25;

phi_min = 0;
phi_max = 11.25;

Eres = 0.22;
E0 = 100;
E1 = E0*(1-Eres);
E2 = E0*(1+Eres);

r_min = sqrt(2*E1);     % inner radius (energy min)
r_max = sqrt(2*E2);     % outer radius (energy max)

% Convert to radians
theta_min = deg2rad(theta_min);
theta_max = deg2rad(theta_max);
phi_min   = deg2rad(phi_min);
phi_max   = deg2rad(phi_max);

%% Resolution (smoothness of surfaces)
n_theta = 40;
n_phi   = 40;
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
surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

% phi_max surface
PHI = phi_max * ones(size(RR));
X = RR .* sin(THETA) .* cos(PHI);
Y = RR .* sin(THETA) .* sin(PHI);
Z = RR .* cos(THETA);
surf(X,Y,Z,'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.8 0.2 0.2]);

%% ---- Plot formatting ----
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3)
camlight; lighting gouraud
grid on
title('Single Spherical Bin')
