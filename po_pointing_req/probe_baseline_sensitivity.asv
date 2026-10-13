%% probe_baseline_sensitivity.m
% Effect of boom geometry on the E-field sensitivity of a spin-plane
% Langmuir-probe pair.
%
% Model: the probes measure dV = -E . L_eff, so the E-noise PSD for a fixed
% voltage-noise PSD scales as
%       S_E = S_V / L_eff^2   ->   S_E(geom) = S_E(ref) * (Lref/Leff)^2
% For the spin-plane component the relevant length is the projection of the
% baseline onto the spin (x-y) plane, Leff = |L_xyz|. The z part of the
% baseline does not help the spin-plane sensitivity; it only matters for
% E_z leakage, shadowing / photoelectron cloud, wake, etc.
%
% Geometry (body frame, x-y = spin plane, z = spin axis):
%   dx    : spacecraft body width along x, ALWAYS present (boom feet at +-dx/2)
%   feet  : f1 = [ dx/2,  dy/2,  dz/2],  f2 = [-dx/2, -dy/2, -dz/2]
%   booms : d1 = [ cos(th), 0, sin(th)],  d2 = [-cos(th), 0, s2*sin(th)]
%           s2 = +1 : both booms angled up ("V")
%           s2 = -1 : one up, one down    ("z-shape", gives dz_tip = 2*Lb*sin(th))
%   probes: p1 = f1 + Lb*d1,  p2 = f2 + Lb*d2,   L = p1 - p2
%
% Case 1: theta = 0, booms along x, feet additionally offset by dy and/or dz
% Case 2: theta > 0, booms angled up/down
% In both cases L_xyz = dx + 2*Lb*cos(th) when dy = 0.

clear; close all;

%% ---- Requirement (spin-plane sensitivity), valid at the reference baseline
f_spec = [1e1 1e2 1e3 1e4 1e5];            % Hz
S_spec = [3e-12 3e-12 2e-12 1e-13 1e-14];  % (V/m)^2/Hz
f      = logspace(1, 5, 300);
S_req  = 10.^interp1(log10(f_spec), log10(S_spec), log10(f), 'linear');

%% ---- Parameters (EDIT)
Lb   = 5;           % boom length [m]
dx   = 1;           % spacecraft body width between boom feet [m]  (placeholder)
Lref = dx + 2*Lb;   % baseline length for which S_spec is assumed [m]
                    % (default: straight booms along x, no y/z offset)

%% ---- Cases to compare
mk = @(name,th,dy,dz,s2) struct('name',name,'dx',dx,'th',th, ...
                                'dy',dy,'dz',dz,'s2',s2);
cases = [ ...
    mk('1a: straight',            0, 0, 0, 1), ...
    mk('1b: dy = 1 m',            0, 1, 0, 1), ...
    mk('1c: dy=1, dz=1 m',        0, 1, 1, 1), ...
    mk('2a: up 20^\circ',        20, 0, 0, 1), ...
    mk('2b: up 30^\circ',        30, 0, 0, 1), ...
    mk('2c: up 45^\circ',        45, 0, 0, 1), ...
    mk('2d: z-shape 30^\circ',   30, 0, 0,-1) ];
N = numel(cases);

%% ---- Evaluate
Lxyz = zeros(1,N); Ltot = zeros(1,N); pen = zeros(1,N);
for k = 1:N
    [p1,p2] = probe_positions(cases(k), Lb);
    L       = p1 - p2;
    Lxyz(k)  = norm(L(1:3));
    Ltot(k) = norm(L); % distance between probes
    pen(k)  = (Lref/Lxyz(k))^2;      % factor on S_E (power)
end
fprintf('dx = %.2f m, Lb = %.2f m, Lref = %.2f m\n', dx, Lb, Lref);
fprintf('%-26s %8s %8s %10s %10s\n','case','Lxyz[m]','|L|[m]','S_E factor','E amp fac');
for k = 1:N
    fprintf('%-26s %8.2f %8.2f %10.2f %10.2f\n', ...
        strrep(cases(k).name,'^\circ','deg'), Lxyz(k), Ltot(k), pen(k), sqrt(pen(k)));
end

%% ---- Fig 1: geometry sketches (top view x-y, side view x-z)
figure('Name','Geometry','Position',[50 50 220*N 450]);
for k = 1:N
    draw_case(cases(k), Lb, k, N, Lxyz(k));
end

%% ---- Fig 2: sensitivity curves vs frequency
figure('Name','Sensitivity vs frequency');
loglog(f, S_req, 'k-', 'LineWidth', 2); hold on
cm = lines(N);
for k = 1:N
    loglog(f, S_req*pen(k), '--', 'Color', cm(k,:), 'LineWidth', 1.3);
end
loglog(f_spec, S_spec, 'ko', 'MarkerFaceColor', 'k');
grid on; xlabel('f [Hz]'); ylabel('S_E [(V/m)^2/Hz]');
title(sprintf('Achieved E-noise PSD for fixed V-noise (spec at L_{ref} = %.1f m)', Lref));
leg = [{'spec @ L_{ref}'}, {cases.name}, {'spec points'}];
legend(leg([1 2:N+1 N+2]), 'Location','southwest', 'Interpreter','tex');

%% ---- Fig 3: penalty map, booms elevated by theta, body width dx
theta_deg = 0:1:80;
dx_vec    = 0.5:0.1:4;
[TH, DX]  = meshgrid(theta_deg, dx_vec);
Lxyz_map   = DX + 2*Lb*cosd(TH);          % same for "V" and z-shape
pen_map   = (Lref ./ Lxyz_map).^2;

figure('Name','Penalty map');
lev = [0.7 0.85 1 1.15 1.5 2 3 5 10];
[C,h] = contourf(TH, DX, pen_map, lev); clabel(C,h); hold on
contour(TH, DX, pen_map, [1 1], 'k', 'LineWidth', 2);
plot(xlim, [dx dx], 'w--');
colorbar; xlabel('boom elevation \theta [deg]'); ylabel('body width d_x [m]');
title('Factor on S_E relative to spec (= (L_{ref}/L_{eff,xyz})^2)');

%% ---- Fig 4: Case 1 - feet offset in y and z, booms along x
dy_vec = 0:0.05:3;
Lxyz_dy = sqrt((dx+2*Lb)^2 + dy_vec.^2);
dz_vec = [0 0.5 1 2];
figure('Name','Case 1');
subplot(1,2,1);
plot(dy_vec, (Lref./Lxyz_dy).^2, 'LineWidth', 1.5); grid on
xlabel('d_y [m]'); ylabel('factor on S_E');
title('Spin-plane penalty (independent of d_z)');
subplot(1,2,2); hold on
for dzv = dz_vec
    plot(dy_vec, abs(dzv)./sqrt((dx+2*Lb)^2 + dy_vec.^2 + dzv^2), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('d_z = %.1f m', dzv));
end
grid on; legend('Location','best'); xlabel('d_y [m]');
ylabel('|L_z| / |L|'); title('Fraction of baseline out of spin plane (E_z leakage)');

%% ---- Fig 5: boom length needed to keep L_xyz = Lref when angled up
th = 0:1:80;
figure('Name','Required boom length'); hold on
for dxv = [1 2 3]
    plot(th, (Lref - dxv)./(2*cosd(th)), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('d_x = %.0f m', dxv));
end
grid on; legend('Location','northwest'); ylim([0 4*Lb]);
xlabel('boom elevation \theta [deg]'); ylabel('L_b needed [m]');
title(sprintf('Boom length for L_{eff,xyz} = %.1f m', Lref));

%% ================= local functions =================
function [p1,p2,f1,f2] = probe_positions(c, Lb)
    f1 = [ c.dx/2,  c.dy/2,  c.dz/2];
    f2 = [-c.dx/2, -c.dy/2, -c.dz/2];
    d1 = [ cosd(c.th), 0,       sind(c.th)];
    d2 = [-cosd(c.th), 0, c.s2*sind(c.th)];
    p1 = f1 + Lb*d1;
    p2 = f2 + Lb*d2;
end

function draw_case(c, Lb, k, N, Lxyz)
    [p1,p2,f1,f2] = probe_positions(c, Lb);
    lim = 1.15*(Lb + c.dx/2);
    subplot(2,N,k); hold on; axis equal; grid on
    plot([f1(1) p1(1)],[f1(2) p1(2)],'b-o');
    plot([f2(1) p2(1)],[f2(2) p2(2)],'r-o');
    plot([f1(1) f2(1)],[f1(2) f2(2)],'k-','LineWidth',3);
    xlim([-lim lim]); ylim([-2 2]);
    xlabel('x [m]'); ylabel('y [m]');
    title({c.name, sprintf('L_{eff,xyz} = %.2f m', Lxyz)}, 'FontSize', 8);
    subplot(2,N,N+k); hold on; axis equal; grid on
    plot([f1(1) p1(1)],[f1(3) p1(3)],'b-o');
    plot([f2(1) p2(1)],[f2(3) p2(3)],'r-o');
    plot([f1(1) f2(1)],[f1(3) f2(3)],'k-','LineWidth',3);
    xlim([-lim lim]); ylim([-Lb Lb]);
    xlabel('x [m]'); ylabel('z [m]');
end
