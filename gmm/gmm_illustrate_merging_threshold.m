% Global normalized L2 merge error for warm + cold two-Gaussian mixtures.
%   Row 1: error maps (3 cold widths) with the err = 0.1 and err = 0.2 curves
%   Row 2: examples on the err = 0.1 curve (one cold width, 4 cold fractions)
%   Row 3: examples on the err = 0.2 curve
%
% Units: warm sigma_w = 1, mean 0. Cold: width sr*sigma_w, offset dmu, fraction r.
% Merge = moment matching (law of total variance); err = ||f-g|| / ||f||.

% ---------------- Settings ----------------
sr_list = [0.1 0.3 0.6];                  % sigma_c / sigma_w
dmu_v   = linspace(0, 4, 61);             % cold offset in units of sigma_w
lr_v    = linspace(-2, log10(0.5), 61);   % log10 of cold fraction r
thrs    = [0.1 0.2];                      % thresholds for the example rows
srEx    = 2;                              % index in sr_list used for the examples
rEx     = [0.2 0.08 0.03];            % cold fractions of the examples

x     = linspace(-8, 14, 6000);           % dx must be << smallest sigma_c
gauss = @(x,m,s) exp(-(x-m).^2./(2*s^2)) ./ (s*sqrt(2*pi));

% ---------------- Scan ----------------
nS = numel(sr_list);  nD = numel(dmu_v);  nR = numel(lr_v);
E  = zeros(nR, nD, nS);
gw = gauss(x, 0, 1);
for is = 1:nS
    s = [1 sr_list(is)];
    for id = 1:nD
        mu = [0 dmu_v(id)];
        gc = gauss(x, mu(2), s(2));
        for ir = 1:nR
            r = 10^lr_v(ir);  w = [1-r r];
            [wM, muM, sM] = merge_gauss(w, mu, s);
            f = (1-r)*gw + r*gc;
            g = wM*gauss(x, muM, sM);
            E(ir,id,is) = sqrt(trapz(x,(f-g).^2) / trapz(x,f.^2));
        end
    end
end

% ---------------- Example points on each threshold curve (chosen cold width) ----------------
nE    = numel(rEx);
rAct  = zeros(1,nE);
dmuEx = nan(numel(thrs), nE);
for k = 1:nE
    [~, ir] = min(abs(lr_v - log10(rEx(k))));
    rAct(k) = 10^lr_v(ir);
    for t = 1:numel(thrs)
        dmuEx(t,k) = cross_dmu(E(ir,:,srEx), dmu_v, thrs(t));
    end
end

% ---------------- Figure ----------------
figure('Position',[50 50 1400 900]);
lineCol = {'r-', 'm--'};   mk = {'o', 's'};

% Row 1: error maps (each panel spans 4 of the 12 columns)
for is = 1:nS
    subplot(3, 3, is); hold on;
    contourf(dmu_v, lr_v, min(E(:,:,is),1), 0:0.05:1, 'LineColor','none');
    for t = 1:numel(thrs)
        contour(dmu_v, lr_v, E(:,:,is), [thrs(t) thrs(t)], lineCol{t}, 'LineWidth', 2);
    end
    if is == srEx
        for t = 1:numel(thrs)
            plot(dmuEx(t,:), log10(rAct), ['k' mk{t}], 'MarkerFaceColor','w', 'MarkerSize', 7);
        end
    end
    caxis([0 1]);
    xlabel('\Delta\mu / \sigma_w');  ylabel('log_{10}(n_{cold}/n_{total})');
    title(sprintf('\\sigma_c/\\sigma_w = %.2f  (red: %.1f, magenta: %.1f)', ...
          sr_list(is), thrs(1), thrs(2)), 'FontSize', 9);
    box on;
    if is == nS, colorbar; end
end

% Rows 2-3: examples (each panel spans 3 of the 12 columns)
fprintf('Examples for sigma_c/sigma_w = %.2f\n', sr_list(srEx));
for t = 1:numel(thrs)
    for k = 1:nE
        subplot(3, 3, 3 + (t-1)*3 + k); hold on;
        if isnan(dmuEx(t,k))
            text(0.5, 0.5, 'no crossing in range', 'HorizontalAlignment','center');
            axis off;  continue
        end
        r  = rAct(k);  s = [1 sr_list(srEx)];  mu = [0 dmuEx(t,k)];  w = [1-r r];
        [wM, muM, sM] = merge_gauss(w, mu, s);
        fw = w(1)*gauss(x, mu(1), s(1));
        fc = w(2)*gauss(x, mu(2), s(2));
        plot(x, fw, '--', 'Color',[0.12 0.47 0.71], 'LineWidth',1.2, 'DisplayName','warm');
        plot(x, fc, '--', 'Color',[0.17 0.63 0.17], 'LineWidth',1.2, 'DisplayName','cold');
        plot(x, fw+fc, 'k-', 'LineWidth',1.8, 'DisplayName','sum (original)');
        plot(x, wM*gauss(x, muM, sM), 'r-', 'LineWidth',1.8, 'DisplayName','merged');
        xlim([-6 10]);  xlabel('x / \sigma_w');  grid on;  box on;
        title(sprintf('thr = %.1f:  r = %.3f, \\Delta\\mu = %.2f\\sigma_w\nT_M/T_c = %.0f', ...
              thrs(t), r, dmuEx(t,k), sM^2/s(2)^2), 'FontSize', 9);
        if k == 1, ylabel('weighted density'); end
        if t == 1 && k == 1, legend('Location','northeast','FontSize',8); end
        fprintf('  thr = %.1f, r = %.3f: dmu = %.2f sigma_w, T_merged/T_cold = %.1f\n', ...
                thrs(t), r, dmuEx(t,k), sM^2/s(2)^2);
    end
end

% ---------------------------------------------------------------------
function d = cross_dmu(row, dmu_v, thr)
    % offset where the error first crosses thr as dmu increases (NaN if it never does)
    j = find((row(1:end-1)-thr).*(row(2:end)-thr) <= 0, 1, 'first');
    if isempty(j), d = NaN; return; end
    den = row(j+1) - row(j);
    if den == 0, t = 0; else, t = (thr - row(j))/den; end
    d = dmu_v(j) + t*(dmu_v(j+1) - dmu_v(j));
end

function [wM, muM, sM] = merge_gauss(w, mu, s)
    wM  = sum(w);
    muM = sum(w.*mu)/wM;
    sM  = sqrt(sum(w.*(s.^2 + (mu-muM).^2))/wM);   % law of total variance
end