% Merging two Gaussians (warm + cold) for different density fractions of the cold one.
% Moment-matched merge (law of total variance); the global normalized L2 error is
% compared with an error normalized by the cold component itself.

mu_w = 0.0;  s_w = 2.5;            % warm, dominant population
mu_c = 1.5;  s_c = 0.2;            % cold population (offset in the same units)
r_cases = [0.4 0.1 0.01];          % cold fraction n_cold / n_total for the example panels
thr = 0.2;                         % example acceptance threshold (your choice)

x     = linspace(-10, 10, 8000);   % fine enough to resolve the cold component
gauss = @(x,m,s) exp(-(x-m).^2./(2*s.^2)) ./ (s*sqrt(2*pi));

% ---------------- Figure 1: three density ratios, linear and log y ----------------
figure('Position',[100 100 1300 650]);
for c = 1:numel(r_cases)
    r  = r_cases(c);
    w  = [1-r, r];  mu = [mu_w mu_c];  s = [s_w s_c];
    [wM, muM, sM, Vw, Vb] = merge_gauss(w, mu, s);

    fw = w(1)*gauss(x,mu(1),s(1));
    fc = w(2)*gauss(x,mu(2),s(2));
    f  = fw + fc;
    g  = wM*gauss(x,muM,sM);
    eG = l2err(x, f, g);                                   % global, as before
    eC = sqrt(trapz(x,(f-g).^2)/trapz(x,fc.^2));           % relative to the cold component

    fprintf('cold fraction %.3f: L2 (global) = %.3f, L2 (vs cold) = %.2f, T_merged/T_cold = %.1f\n', ...
            r, eG, eC, sM^2/s_c^2);

    for row = 1:2
        subplot(2,3,(row-1)*3+c); hold on;
        plot(x, fw, '--', 'Color',[0.12 0.47 0.71], 'LineWidth',1.2, 'DisplayName','warm');
        plot(x, fc, '--', 'Color',[0.17 0.63 0.17], 'LineWidth',1.2, 'DisplayName','cold');
        plot(x, f,  'k-', 'LineWidth',1.8, 'DisplayName','sum (original)');
        plot(x, g,  'r-', 'LineWidth',1.8, 'DisplayName','merged');
        xlabel('x'); ylabel('weighted density'); grid on; box on;
        if row == 1
            title(sprintf('cold fraction = %.2f\nL2 global = %.3f,  L2 vs cold = %.2f', r, eG, eC));
            legend('Location','northeast');
        else
            set(gca,'YScale','log'); ylim([1e-6 1]);
            title('same, log y');
        end
    end
end

% ---------------- Figure 2: scan the cold fraction ----------------
rs = logspace(-3, log10(0.5), 40);
eG = zeros(size(rs));  eC = eG;  Tr = eG;
for n = 1:numel(rs)
    r = rs(n);  w = [1-r, r];  mu = [mu_w mu_c];  s = [s_w s_c];
    [wM, muM, sM] = merge_gauss(w, mu, s);
    fc = w(2)*gauss(x,mu(2),s(2));
    f  = w(1)*gauss(x,mu(1),s(1)) + fc;
    g  = wM*gauss(x,muM,sM);
    eG(n) = l2err(x, f, g);
    eC(n) = sqrt(trapz(x,(f-g).^2)/trapz(x,fc.^2));
    Tr(n) = sM^2/s_c^2;
end

figure('Position',[100 100 1000 400]);
subplot(1,2,1); hold on;
plot(rs, eG, 'k-',  'LineWidth',1.5, 'DisplayName','normalized by ||f|| (global)');
plot(rs, eC, 'r-',  'LineWidth',1.5, 'DisplayName','normalized by ||cold component||');
yline(thr, 'k--', 'DisplayName','example threshold');
set(gca,'XScale','log','YScale','log');
xlabel('cold fraction  n_{cold}/n_{total}'); ylabel('L2 error');
title('Merge error vs. density of the cold component');
legend('Location','best'); grid on; box on;

subplot(1,2,2);
loglog(rs, Tr, 'b-', 'LineWidth',1.5); grid on; box on;
xlabel('cold fraction  n_{cold}/n_{total}');
ylabel('T_{merged} / T_{cold}');
title('Temperature of the merged component');

% ---------------------------------------------------------------------
function [wM, muM, sM, Vw, Vb] = merge_gauss(w, mu, s)
    wM  = sum(w);
    muM = sum(w.*mu)/wM;
    Vw  = sum(w.*s.^2)/wM;
    Vb  = sum(w.*(mu-muM).^2)/wM;
    sM  = sqrt(Vw + Vb);
end

function e = l2err(x, f, g)
    e = sqrt(trapz(x,(f-g).^2) / trapz(x,f.^2));
end