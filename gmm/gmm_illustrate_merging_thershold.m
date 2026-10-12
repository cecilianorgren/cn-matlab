% Global normalized L2 merge error, warm + cold Gaussian, for sigma_c/sigma_w = 0.1.
%   Top row:    error map (offset vs. cold fraction) with the err = 0.1 and 0.2 curves
%   Rows 2, 3:  examples at given offsets dmu; for each offset the cold fraction r is
%               solved so that err = threshold (0.1 in row 2, 0.2 in row 3)
%
% Units: warm sigma_w = 1, mean 0. Cold: width sr*sigma_w, offset dmu, fraction r.
% Merge = moment matching (law of total variance); err = ||f-g|| / ||f||.

% ---------------- Settings ----------------
sr     = 0.3;                             % sigma_c / sigma_w
dmu_ex = [1 2];                       % offsets of the examples (units of sigma_w)
thrs   = [0.1 0.2];
dmu_v  = linspace(0, 3, 61);              % map axes
lr_v   = linspace(-2, log10(0.5), 61);    % log10 of cold fraction
lr_scan = linspace(-3, log10(0.5), 200);  % scan used to bracket the crossing in r

x     = linspace(-8, 14, 6000);           % dx must be << sr
gauss = @(x,m,s) exp(-(x-m).^2./(2*s^2)) ./ (s*sqrt(2*pi));
errfun = @(lr, dmu) err_two(x, gauss, lr, dmu, sr);   % error for given log10(r) and offset

% ---------------- Error map ----------------
nD = numel(dmu_v);  nR = numel(lr_v);
E  = zeros(nR, nD);
for id = 1:nD
    for ir = 1:nR
        E(ir,id) = errfun(lr_v(ir), dmu_v(id));
    end
end

% ---------------- Cold fraction at err = thr, for each example offset ----------------
nE = numel(dmu_ex);
lrEx = nan(numel(thrs), nE);
for t = 1:numel(thrs)
    for k = 1:nE
        lrEx(t,k) = find_lr(dmu_ex(k), thrs(t), lr_scan, errfun);
    end
end

% ---------------- Figure ----------------
%figure('Position',[50 50 1300 850]);
figure
lineCol = {'r-', 'm--'};  mk = {'o', 's'};

subplot(3, 2, 1:2); hold on;
contourf(dmu_v, lr_v, min(E,1), 0:0.05:1, 'LineColor','none');
for t = 1:numel(thrs)
    contour(dmu_v, lr_v, E, [thrs(t) thrs(t)], lineCol{t}, 'LineWidth', 2);
    plot(dmu_ex, lrEx(t,:), ['k' mk{t}], 'MarkerFaceColor','w', 'MarkerSize', 7);
end
caxis([0 1]);  colorbar;col
xlabel('\Delta\mu / \sigma_w');  ylabel('log_{10}(n_{cold}/n_{total})');
title(sprintf('\\sigma_c/\\sigma_w = %.2f   (red: err = %.1f, magenta: err = %.1f)', sr, thrs(1), thrs(2)));
box on;

fprintf('sigma_c/sigma_w = %.2f\n', sr);
for t = 1:numel(thrs)
    for k = 1:nE
        subplot(3, 2, 2*t + k); hold on;
        if isnan(lrEx(t,k))
            text(0.5, 0.5, 'no crossing in range', 'HorizontalAlignment','center');
            axis off;  continue
        end
        r  = 10^lrEx(t,k);  s = [1 sr];  mu = [0 dmu_ex(k)];  w = [1-r r];
        [wM, muM, sM] = merge_gauss(w, mu, s);
        fw = w(1)*gauss(x, mu(1), s(1));
        fc = w(2)*gauss(x, mu(2), s(2));
        plot(x, fw, '--', 'Color',[0.12 0.47 0.71], 'LineWidth',1.2, 'DisplayName','warm');
        plot(x, fc, '--', 'Color',[0.17 0.63 0.17], 'LineWidth',1.2, 'DisplayName','cold');
        plot(x, fw+fc, 'k-', 'LineWidth',1.8, 'DisplayName','sum (original)');
        plot(x, wM*gauss(x, muM, sM), 'r-', 'LineWidth',1.8, 'DisplayName','merged');
        xlim([-6 8]);  xlabel('x / \sigma_w');  grid on;  box on;
        title(sprintf('err = %.1f:  \\Delta\\mu = %.1f\\sigma_w, r = %.4f\nT_M/T_c = %.0f', ...
              thrs(t), dmu_ex(k), r, sM^2/s(2)^2), 'FontSize', 9);
        if k == 1, ylabel('weighted density'); end
        if t == 1 && k == 1, legend('Location','northeast','FontSize',8,'Box','off'); end
        fprintf('  err = %.1f, dmu = %.1f: r = %.4f, T_merged/T_cold = %.1f\n', ...
                thrs(t), dmu_ex(k), r, sM^2/s(2)^2);
    end
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
h(1).Position(1) = 0.3;
h(1).Position(3) = 0.4;
c_eval('h(?).XLim = [-4 4];',2:5)
% ---------------------------------------------------------------------
function e = err_two(x, gauss, lr, dmu, sr)
    r = 10^lr;  w = [1-r r];  mu = [0 dmu];  s = [1 sr];
    [wM, muM, sM] = merge_gauss(w, mu, s);
    f = w(1)*gauss(x,mu(1),s(1)) + w(2)*gauss(x,mu(2),s(2));
    g = wM*gauss(x, muM, sM);
    e = sqrt(trapz(x,(f-g).^2) / trapz(x,f.^2));
end

function lx = find_lr(dmu, thr, lr_scan, errfun)
    % log10(r) at which err first reaches thr as r increases (NaN if it never does)
    e = arrayfun(@(l) errfun(l, dmu), lr_scan);
    j = find(e(1:end-1) < thr & e(2:end) >= thr, 1, 'first');
    if isempty(j), lx = NaN; return; end
    lx = fzero(@(l) errfun(l, dmu) - thr, [lr_scan(j) lr_scan(j+1)]);
end

function [wM, muM, sM] = merge_gauss(w, mu, s)
    wM  = sum(w);
    muM = sum(w.*mu)/wM;
    sM  = sqrt(sum(w.*(s.^2 + (mu-muM).^2))/wM);   % law of total variance
end