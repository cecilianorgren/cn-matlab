% Law of total variance for a two-component Gaussian mixture (1D).
% Merge two components into one by matching total weight, mean and variance,
% then measure how well the merged Gaussian reproduces the original sum
% with a normalized L2 error.

% --- Two components ---
w  = [0.6 0.4];          % weights
mu = [0.0 1.5];          % means
s  = [1.0 0.6];          % standard deviations

gauss = @(x,m,sg) exp(-(x-m).^2./(2*sg.^2)) ./ (sg*sqrt(2*pi));

% --- Merge (moment matching) ---
merge = @(w,mu,s) deal( sum(w), ...
                        sum(w.*mu)/sum(w), ...
                        sqrt( sum(w.*(s.^2 + (mu - sum(w.*mu)/sum(w)).^2)) / sum(w) ) );
[wM, muM, sM] = merge(w,mu,s);

% Terms of the law of total variance
mu_tot   = sum(w.*mu)/sum(w);
V_within = sum(w.*s.^2)/sum(w);
V_betw   = sum(w.*(mu-mu_tot).^2)/sum(w);
fprintf('within  = %.4f\nbetween = %.4f\ntotal   = %.4f  (merged sigma^2 = %.4f)\n', ...
        V_within, V_betw, V_within+V_betw, sM^2);

% --- Normalized L2 error: ||f - g|| / ||f||, f = original sum, g = merged ---
x     = linspace(-8, 10, 4000);
l2err = @(f,g) sqrt(trapz(x,(f-g).^2) / trapz(x,f.^2));

f1 = w(1)*gauss(x,mu(1),s(1));
f2 = w(2)*gauss(x,mu(2),s(2));
f  = f1 + f2;
g  = wM*gauss(x,muM,sM);
err = l2err(f,g);

% Check that the sum and the merged Gaussian really share mean and variance
m_f = trapz(x,x.*f)/trapz(x,f);
v_f = trapz(x,(x-m_f).^2.*f)/trapz(x,f);
fprintf('sum:    mean = %.4f, var = %.4f\nmerged: mean = %.4f, var = %.4f\n', m_f, v_f, muM, sM^2);
fprintf('normalized L2 error = %.3f\n', err);

% --- Scan the separation of the two means to see how the error behaves ---
sep   = linspace(0, 6, 61);             % (mu2 - mu1) / mean sigma
sbar  = mean(s);
errS  = zeros(size(sep));
for n = 1:numel(sep)
    mu_n = [0, sep(n)*sbar];
    [wm,mm,sm] = merge(w,mu_n,s);
    fn = w(1)*gauss(x,mu_n(1),s(1)) + w(2)*gauss(x,mu_n(2),s(2));
    gn = wm*gauss(x,mm,sm);
    errS(n) = l2err(fn,gn);
end
thr = 0.1;                              % example acceptance threshold (your choice)

% --- Plot ---
figure('Position',[100 100 1000 400]);

subplot(1,2,1); hold on;
plot(x, f1, '--', 'Color', [0.12 0.47 0.71], 'LineWidth', 1.2, 'DisplayName','component 1');
plot(x, f2, '--', 'Color', [0.17 0.63 0.17], 'LineWidth', 1.2, 'DisplayName','component 2');
plot(x, f,  'k-',  'LineWidth', 1.8, 'DisplayName','sum (original)');
plot(x, g,  'r-',  'LineWidth', 1.8, 'DisplayName','merged Gaussian');
xlabel('x'); ylabel('weighted density');
title(sprintf('Merge: normalized L2 error = %.3f', err));
legend('Location','northwest','Box','off'); grid on; box on;
xlim([-5 5])

subplot(1,2,2); hold on;
plot(sep, errS, 'k-', 'LineWidth', 1.5, 'DisplayName','L2 error');
yline(thr, 'r--', 'DisplayName','example threshold');
plot((mu(2)-mu(1))/sbar, err, 'ro', 'MarkerFaceColor','r', 'DisplayName','current case');
xlabel('(\mu_2-\mu_1) / mean \sigma'); ylabel('normalized L2 error');
title('Error vs. separation of the means');
legend('Location','northwest','Box','off'); grid on; box on;

%%
