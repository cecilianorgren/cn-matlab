% Law of total variance: merge two Gaussians by moment matching, for a
% "good" merge (overlapping components) and a "bad" merge (separated components).
% Both cases share weights and widths, so they differ only in the separation of the means.

w  = [0.6 0.4];
s  = [1.0 0.6];
cases(1).name = 'good merge (overlapping)'; 
cases(1).mu = [0.0 0.8];
cases(1).w = [0.4 0.6];
cases(2).name = 'bad merge (separated)';    
cases(2).mu = [0.0 2.5];
cases(1).w = [0.2 0.8];

x   = linspace(-8, 10, 4000);
thr = 0.2;                                  % example acceptance threshold (your choice)

gauss = @(x,m,sg) exp(-(x-m).^2./(2*sg.^2)) ./ (sg*sqrt(2*pi));

%figure('Position',[100 100 1300 400]);
figure(43)
col = [0.12 0.47 0.71; 0.17 0.63 0.17];
sbar = mean(s);

for c = 1:numel(cases)
    mu = cases(c).mu;
    [wM, muM, sM, Vw, Vb] = merge_gauss(w, mu, s);

    f1 = w(1)*gauss(x,mu(1),s(1));
    f2 = w(2)*gauss(x,mu(2),s(2));
    f  = f1 + f2;
    g  = wM*gauss(x,muM,sM);
    cases(c).err = l2err(x,f,g);

    fprintf('%s:\n  within = %.3f, between = %.3f, total = %.3f (merged sigma^2 = %.3f)\n', ...
            cases(c).name, Vw, Vb, Vw+Vb, sM^2);
    fprintf('  normalized L2 error = %.3f  -> %s\n', cases(c).err, ...
            ternary(cases(c).err < thr, 'accept', 'reject'));

    subplot(1,3,c); hold on;
    plot(x, f1, '--', 'Color', col(1,:), 'LineWidth', 1.2, 'DisplayName','comp. 1');
    plot(x, f2, '--', 'Color', col(2,:), 'LineWidth', 1.2, 'DisplayName','comp. 2');
    plot(x, f,  'k-', 'LineWidth', 1.8, 'DisplayName','sum (original)');
    plot(x, g,  'r-', 'LineWidth', 1.8, 'DisplayName','merged Gauss.');
    xlabel('x'); ylabel('weighted density');
    title(sprintf('%s\nL2 error = %.3f', cases(c).name, cases(c).err));
    legend('Location','northwest','Box','off'); grid on; box on;
    xlim([-7 7])
end

% --- Error vs. separation of the means (same w and s as above) ---
sep  = linspace(0, 5, 61);
errS = zeros(size(sep));
for n = 1:numel(sep)
    mu_n = [0, sep(n)*sbar];
    [wm,mm,sm] = merge_gauss(w, mu_n, s);
    fn = w(1)*gauss(x,mu_n(1),s(1)) + w(2)*gauss(x,mu_n(2),s(2));
    errS(n) = l2err(x, fn, wm*gauss(x,mm,sm));
end

% --- Error vs. separation of the means (same w and s as above) ---
weightratio  = linspace(0.1, 5, 61);
errSw = zeros(size(weightratio));
for n = 1:numel(weightratio)
    w_n = w(1)*[1, weightratio(n)];
    [wm,mm,sm] = merge_gauss(w_n, mu, s);
    fn = w_n(1)*gauss(x,mu(1),s(1)) + w_n(2)*gauss(x,mu(2),s(2));
    errSw(n) = l2err(x, fn, wm*gauss(x,mm,sm));
end

subplot(1,3,3); hold on;

plot(weightratio, errSw, 'k-', 'LineWidth', 1.5, 'DisplayName','L2 error');
yline(thr, 'r--', 'DisplayName','example threshold');
for c = 1:numel(cases)
    w_p = w(c)/w(1);%(cases(c).w(2)-cases(c).w(1));
    plot(w_p, cases(c).err, 'o', 'MarkerFaceColor', [c==2 0.5*(c==1) 0], ...
         'MarkerEdgeColor','k', 'DisplayName', cases(c).name);
end
xlabel('(\mu_2-\mu_1) / mean \sigma'); ylabel('normalized L2 error');
title('Error vs. separation of the means');
legend('Location','northwest','Box','off'); grid on; box on;

%%
plot(sep, errS, 'k-', 'LineWidth', 1.5, 'DisplayName','L2 error');
yline(thr, 'r--', 'DisplayName','example threshold');
for c = 1:numel(cases)
    d = (cases(c).mu(2)-cases(c).mu(1))/sbar;
    plot(d, cases(c).err, 'o', 'MarkerFaceColor', [c==2 0.5*(c==1) 0], ...
         'MarkerEdgeColor','k', 'DisplayName', cases(c).name);
end
xlabel('(\mu_2-\mu_1) / mean \sigma'); ylabel('normalized L2 error');
title('Error vs. separation of the means');
legend('Location','northwest','Box','off'); grid on; box on;

% ---------------------------------------------------------------------
function [wM, muM, sM, Vw, Vb] = merge_gauss(w, mu, s)
    wM  = sum(w);
    muM = sum(w.*mu)/wM;
    Vw  = sum(w.*s.^2)/wM;               % average within-component variance
    Vb  = sum(w.*(mu-muM).^2)/wM;        % between-component spread of the means
    sM  = sqrt(Vw + Vb);                 % law of total variance
end

function e = l2err(x, f, g)
    e = sqrt(trapz(x,(f-g).^2) / trapz(x,f.^2));   % ||f-g|| / ||f||
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end