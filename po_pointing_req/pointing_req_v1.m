%% probe_baseline_sensitivity.m
% Effect of boom geometry on the E-field sensitivity of a spin-plane
% Langmuir-probe pair.
%
% Model: the probes measure dV = -E . L_eff, so the E-noise PSD for a fixed
% voltage-noise PSD scales as
%       S_E = S_V / L_eff^2   ->   S_E(geom) = S_E(ref) * (Lref/Leff)^2
% For the spin-plane component the relevant length is the projection of the
% baseline onto the spin (x-y) plane, Leff = |L_xy|. The z part of the
% baseline does not help the spin-plane sensitivity; it only matters for
% E_z leakage, shadowing / photoelectron cloud, wake, etc.
%
% Geometry (body frame, x-y = spin plane, z = spin axis):
%   feet  : f1 = [ Lsc/2,  dy/2,  dz/2],  f2 = [-Lsc/2, -dy/2, -dz/2]
%   booms : d1 = [ cos(th), 0, sin(th)],  d2 = [-cos(th), 0, s2*sin(th)]
%           s2 = +1 : both booms angled up ("V")
%           s2 = -1 : one up, one down    ("z-shape", gives dz = 2*Lb*sin(th))
%   probes: p1 = f1 + Lb*d1,  p2 = f2 + Lb*d2,   L = p1 - p2
%
% Case 1: theta = 0, feet separated by dy and/or dz (Lsc = body width along x)
% Case 2: theta > 0, feet separated by Lsc along x

clear; close all;

%% ---- Requirement (spin-plane sensitivity), valid at the reference baseline
f_spec = [1e1 1e2 1e3 1e4 1e5];            % Hz
S_spec = [3e-12 3e-12 2e-12 1e-13 1e-14];  % (V/m)^2/Hz
f      = logspace(1, 5, 300);
S_req  = 10.^interp1(log10(f_spec), log10(S_spec), log10(f), 'linear');

%% ---- Parameters (EDIT)
Lsc  = 1;
Lb   = 5;        % boom length [m]
Lref = 2*Lb;     % baseline length for which S_spec is assumed [m]
                 % (default: straight, no foot separation)

%% ---- Cases to compare
mk = @(name,Lsc,th,dy,dz,s2) struct('name',name,'Lsc',Lsc,'th',th, ...
                                    'dy',dy,'dz',dz,'s2',s2);
cases = [ ...
    mk('1a: straight',            0,   0, 0, 0, 1), ...
    mk('1b: dy = 1 m',            0,   0, 1, 0, 1), ...
    mk('1c: dy=1, dz=1 m',        0,   0, 1, 1, 1), ...
    mk('2a: up 30^\circ, Lsc=0',  0,  30, 0, 0, 1), ...
    mk('2b: up 30^\circ, Lsc=1',  1,  30, 0, 0, 1), ...
    mk('2c: up 45^\circ, Lsc=1',  1,  45, 0, 0, 1), ...
    mk('2d: z-shape 30^\circ',    1,  30, 0, 0,-1) ];
N = numel(cases);

%% ---- Evaluate
Lxy = zeros(1,N); Ltot = zeros(1,N); pen = zeros(1,N);
for k = 1:N
    [p1,p2] = probe_positions(cases(k), Lb);
    L       = p1 - p2;
    Lxy(k)  = norm(L(1:2));
    Ltot(k) = norm(L);
    pen(k)  = (Lref/Lxy(k))^2;      % factor on S_E (power)
end
fprintf('%-26s %8s %8s %10s %10s\n','case','Lxy[m]','|L|[m]','S_E factor','E amp fac');
for k = 1:N
    fprintf('%-26s %8.2f %8.2f %10.2f %10.2f\n', ...
        strrep(cases(k).name,'^\circ','deg'), Lxy(k), Ltot(k), pen(k), sqrt(pen(k)));
end

%% ---- Fig 1: geometry sketches (top view x-y, side view x-z)
figure('Name','Geometry','Position',[50 50 220*N 450]);
for k = 1:N
    draw_case(cases(k), Lb, k, N, Lxy(k));
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
leg = [{'spec @ L_{ref}'}, strrep({cases.name},'_','\_')];
legend(leg, 'Location','southwest', 'Interpreter','tex');

%% ---- Fig 3: penalty map, both booms up by theta, foot separation Lsc
theta_deg = 0:1:80;
Lsc_vec   = 0:0.1:3;
[TH, LSC] = meshgrid(theta_deg, Lsc_vec);
Lxy_map   = LSC + 2*Lb*cosd(TH);          % same for "V" and z-shape
pen_map   = (Lref ./ Lxy_map).^2;

figure('Name','Penalty map');
lev = [0.5 0.7 0.85 1 1.15 1.5 2 3 5 10];
[C,h] = contourf(TH, LSC, pen_map, lev); clabel(C,h); hold on
contour(TH, LSC, pen_map, [1 1], 'k', 'LineWidth', 2);
colorbar; xlabel('boom elevation \theta [deg]'); ylabel('L_{sc} (foot separation) [m]');
title('Factor on S_E relative to spec (= (L_{ref}/L_{eff,xy})^2)');

%% ---- Fig 4: Case 1 - feet offset in y and z, booms along x
dy_vec = 0:0.05:3;
Lxy_dy = sqrt((2*Lb)^2 + dy_vec.^2);
dz_vec = [0 0.5 1 2];
figure('Name','Case 1');
subplot(1,2,1);
plot(dy_vec, (Lref./Lxy_dy).^2, 'LineWidth', 1.5); grid on
xlabel('d_y [m]'); ylabel('factor on S_E');
title('Spin-plane penalty (independent of d_z)');
subplot(1,2,2); hold on
for dz = dz_vec
    plot(dy_vec, abs(dz)./sqrt((2*Lb)^2 + dy_vec.^2 + dz^2), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('d_z = %.1f m', dz));
end
grid on; legend('Location','best'); xlabel('d_y [m]');
ylabel('|L_z| / |L|'); title('Fraction of baseline out of spin plane (E_z leakage)');

%% ================= local functions =================
function [p1,p2,f1,f2] = probe_positions(c, Lb)
    f1 = [ c.Lsc/2,  c.dy/2,  c.dz/2];
    f2 = [-c.Lsc/2, -c.dy/2, -c.dz/2];
    d1 = [ cosd(c.th), 0,       sind(c.th)];
    d2 = [-cosd(c.th), 0, c.s2*sind(c.th)];
    p1 = f1 + Lb*d1;
    p2 = f2 + Lb*d2;
end

function draw_case(c, Lb, k, N, Lxy)
    [p1,p2,f1,f2] = probe_positions(c, Lb);
    lim = 1.15*(Lb + c.Lsc/2);
    subplot(2,N,k); hold on; axis equal; grid on
    plot([f1(1) p1(1)],[f1(2) p1(2)],'b-o');
    plot([f2(1) p2(1)],[f2(2) p2(2)],'r-o');
    plot([f1(1) f2(1)],[f1(2) f2(2)],'k-','LineWidth',3);
    xlim([-lim lim]); ylim([-2 2]);
    xlabel('x [m]'); ylabel('y [m]');
    title({c.name, sprintf('L_{eff,xy} = %.2f m', Lxy)}, 'FontSize', 8);
    subplot(2,N,N+k); hold on; axis equal; grid on
    plot([f1(1) p1(1)],[f1(3) p1(3)],'b-o');
    plot([f2(1) p2(1)],[f2(3) p2(3)],'r-o');
    plot([f1(1) f2(1)],[f1(3) f2(3)],'k-','LineWidth',3);
    xlim([-lim lim]); ylim([-Lb Lb]);
    xlabel('x [m]'); ylabel('z [m]');
end