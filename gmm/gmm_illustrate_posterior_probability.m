rng(42);

N = 1e5;                          % same number of particles in both populations
sigma_cold = 0.5;
sigma_warm = 2.0;

cold = sigma_cold * randn(N, 1);
warm = sigma_warm * randn(N, 1);

edges = linspace(-8, 8, 161);

figure; hold on;

% 'pdf' normalization -> area under each histogram is 1
histogram(cold, edges, 'Normalization', 'pdf', 'FaceColor', [0.12 0.47 0.71], ...
  'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
  'DisplayName', sprintf('cold (\\sigma = %.1f)', sigma_cold));
histogram(warm, edges, 'Normalization', 'pdf', 'FaceColor', [0.84 0.15 0.16], ...
  'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
  'DisplayName', sprintf('warm (\\sigma = %.1f)', sigma_warm));

% analytic Gaussians for reference
x = linspace(-8, 8, 800);
gauss = @(x, s) exp(-x.^2 / (2*s^2)) / (s * sqrt(2*pi));
plot(x, gauss(x, sigma_cold), 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(x, gauss(x, sigma_warm), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('x (or v)');
ylabel('probability density');
title('Same N, different variance: normalization forces a taller peak for cold');
legend('show');
box on;
hold off;

%%
rng(42);

N = 1e5;                          % same number of particles in both populations
sigma_cold = 0.5;
sigma_warm = 2.0;

cold = sigma_cold * randn(N, 1);
warm = sigma_warm * randn(N, 1);

gauss = @(x, m, s) exp(-(x-m).^2 / (2*s.^2)) ./ (s * sqrt(2*pi));

% --- Fit a Gaussian to each population (MLE) and compute logL and BIC ---
k = 2;                            % free parameters: mean, sigma
fitStats = @(d) deal(mean(d), std(d, 1));   % MLE mean and sigma (normalize by N)

[mu_c, s_c] = fitStats(cold);
[mu_w, s_w] = fitStats(warm);

logL_cold = sum(log(gauss(cold, mu_c, s_c)));
logL_warm = sum(log(gauss(warm, mu_w, s_w)));

BIC_cold = k*log(N) - 2*logL_cold;
BIC_warm = k*log(N) - 2*logL_warm;

fprintf('          logL          BIC\n');
fprintf('cold: %12.1f %12.1f\n', logL_cold, BIC_cold);
fprintf('warm: %12.1f %12.1f\n', logL_warm, BIC_warm);

% Analytic expectation: logL/N = -0.5*log(2*pi*e*sigma^2)  (minus the entropy)
fprintf('\nlogL/N  cold: %.3f (analytic %.3f)\n', logL_cold/N, -0.5*log(2*pi*exp(1)*sigma_cold^2));
fprintf('logL/N  warm: %.3f (analytic %.3f)\n', logL_warm/N, -0.5*log(2*pi*exp(1)*sigma_warm^2));

% --- Plot ---
edges = linspace(-8, 8, 161);
figure; hold on;

histogram(cold, edges, 'Normalization', 'pdf', 'FaceColor', [0.12 0.47 0.71], ...
  'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
  'DisplayName', sprintf('cold (\\sigma = %.1f)', sigma_cold));
histogram(warm, edges, 'Normalization', 'pdf', 'FaceColor', [0.84 0.15 0.16], ...
  'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
  'DisplayName', sprintf('warm (\\sigma = %.1f)', sigma_warm));

x = linspace(-8, 8, 800);
plot(x, gauss(x, mu_c, s_c), 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(x, gauss(x, mu_w, s_w), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('x (or v)');
ylabel('probability density');
title('Both fits are perfect, yet logL and BIC differ purely because of the width');
legend('show');

txt = sprintf(['cold: logL = %.0f,  BIC = %.0f\n' ...
  'warm: logL = %.0f,  BIC = %.0f'], ...
  logL_cold, BIC_cold, logL_warm, BIC_warm);
text(0.98, 0.70, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
  'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

box on; hold off;