% Merge error of a two-Gaussian mixture (warm + cold beam) over a parameter space,
% and example mixtures along the err = threshold curve.
%
% Warm component: sigma_w = 1, mean 0 (everything is in units of sigma_w).
% Cold component: mean offset dmu (in sigma_w), width ratio sr = sigma_c/sigma_w,
%                 density fraction r = n_c/n_total.
% Merge = moment matching (law of total variance).
% Two error measures:
%   1) global: ||f-g|| / ||f||          (f = original sum, g = merged Gaussian)
%   2) cold:   ||f-g|| / ||n_c f_c||    (is the cold population still represented?)

% ---------------- Settings ----------------
sr_list   = [0.1 0.3 0.6];             % sigma_c / sigma_w  (T_c/T_w = sr^2)
dmu_v     = linspace(0, 6, 61);        % offset of the cold mean, units of sigma_w
lr_v      = linspace(-3, log10(0.5), 61);   % log10 of cold fraction r
thr       = 0.1;                       % acceptance threshold (your choice)
metricSel = 1;                         % metric used for the examples: 1 = global, 2 = cold
srEx      = 2;                         % index in sr_list used for the example figure
rEx       = [0.3 0.1 0.03 0.01];       % cold fractions at which examples are taken

x     = linspace(-8, 14, 6000);        % dx must be << smallest sigma (here 0.1)
gauss = @(x,m,s) exp(-(x-m).^2./(2*s^2)) ./ (s*sqrt(2*pi));

% ---------------- Scan ----------------
nS = numel(sr_list);  nD = numel(dmu_v);  nR = numel(lr_v);
E  = zeros(nR, nD, nS, 2);
gw = gauss(x, 0, 1);
for is = 1:nS
    s = [1 sr_list(is)];
    for id = 1:nD
        mu = [0 dmu_v(id)];
        gc = gauss(x, mu(2), s(2));
        for ir = 1:nR
            r  = 10^lr_v(ir);
            w  = [1-r r];
            [wM, muM, sM] = merge_gauss(w, mu, s);
            fc = r*gc;
            f  = (1-r)*gw + fc;
            g  = wM*gauss(x, muM, sM);
            d2 = trapz(x, (f-g).^2);
            E(ir,id,is,1) = sqrt(d2/trapz(x, f.^2));
            E(ir,id,is,2) = sqrt(d2/trapz(x, fc.^2));
        end
    end
end

% ---------------- Points on the err = thr curve (for the chosen sigma ratio) ----------------
Em    = E(:,:,srEx,metricSel);
rAct  = nan(size(rEx));  dmuEx = nan(size(rEx));
for k = 1:numel(rEx)
    [~, ir] = min(abs(lr_v - log10(rEx(k))));
    rAct(k) = 10^lr_v(ir);
    row = Em(ir,:);
    j = find((row(1:end-1)-thr).*(row(2:end)-thr) <= 0, 1, 'first');   % first crossing in dmu
    if ~isempty(j)
        den = row(j+1) - row(j);
        t   = (den ~= 0) * (thr-row(j)) / (den + (den==0));
        dmuEx(k) = dmu_v(j) + t*(dmu_v(j+1)-dmu_v(j));
    end
end

% ---------------- Figure 1: error maps ----------------
names  = {'global  ||f-g||/||f||', 'cold  ||f-g||/||n_c f_c||'};
levels = 0:0.05:1;
figure('Position',[50 50 1300 700]);
for m = 1:2
    for is = 1:nS
        subplot(2, nS, (m-1)*nS + is);
        contourf(dmu_v, lr_v, min(E(:,:,is,m),1), levels, 'LineColor','none'); hold on;
        contour(dmu_v, lr_v, E(:,:,is,m), [thr thr], 'r-', 'LineWidth', 2);
        if m == metricSel && is == srEx
            plot(dmuEx, log10(rAct), 'wo', 'MarkerFaceColor','k', 'MarkerSize', 7);
        end
        clim_ = [0 1]; caxis(clim_);
        xlabel('\Delta\mu / \sigma_w');  ylabel('log_{10}(n_{cold}/n_{total})');
        title(sprintf('%s\n\\sigma_c/\\sigma_w = %.2f,  red: err = %.2f', names{m}, sr_list(is), thr));
        box on;
        if is == nS, colorbar; end
    end
end

% ---------------- Figure 2: examples along the threshold curve ----------------
nE = numel(rEx);
figure('Position',[50 50 1300 600]);
fprintf('Examples on err = %.2f (%s metric), sigma_c/sigma_w = %.2f\n', thr, names{metricSel}, sr_list(srEx));
for k = 1:nE
    for row = 1:2
        subplot(2, nE, (row-1)*nE + k); hold on;
        if isnan(dmuEx(k))
            text(0.5, 0.5, 'no crossing in range', 'HorizontalAlignment','center');
            axis off;  continue
        end
        r  = rAct(k);  s = [1 sr_list(srEx)];  mu = [0 dmuEx(k)];  w = [1-r r];
        [wM, muM, sM] = merge_gauss(w, mu, s);
        fw = w(1)*gauss(x,mu(1),s(1));
        fc = w(2)*gauss(x,mu(2),s(2));
        plot(x, fw,        '--', 'Color',[0.12 0.47 0.71], 'LineWidth',1.2, 'DisplayName','warm');
        plot(x, fc,        '--', 'Color',[0.17 0.63 0.17], 'LineWidth',1.2, 'DisplayName','cold');
        plot(x, fw+fc,     'k-', 'LineWidth',1.8, 'DisplayName','sum');
        plot(x, wM*gauss(x,muM,sM), 'r-', 'LineWidth',1.8, 'DisplayName','merged');
        xlabel('x / \sigma_w'); grid on; box on; xlim([-6 10]);
        if row == 1
            title(sprintf('r = %.3f, \\Delta\\mu = %.2f\\sigma_w\nT_M/T_c = %.0f', r, dmuEx(k), sM^2/s(2)^2));
            ylabel('weighted density');
            if k == 1, legend('Location','northeast'); end
        else
            set(gca,'YScale','log'); ylim([1e-6 3]); ylabel('same, log y');
        end
        if row == 1
            fprintf('  r = %.3f: dmu = %.2f sigma_w, T_merged/T_cold = %.1f\n', r, dmuEx(k), sM^2/s(2)^2);
        end
    end
end

% ---------------------------------------------------------------------
function [wM, muM, sM, Vw, Vb] = merge_gauss(w, mu, s)
    wM  = sum(w);
    muM = sum(w.*mu)/wM;
    Vw  = sum(w.*s.^2)/wM;               % within-component variance
    Vb  = sum(w.*(mu-muM).^2)/wM;        % between-component spread of the means
    sM  = sqrt(Vw + Vb);                 % law of total variance
end