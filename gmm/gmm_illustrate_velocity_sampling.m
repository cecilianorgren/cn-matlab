% Effect of the sampling variable (uniform in v^p) on the moments of a
% particle population drawn from one energy bin, vs. energy resolution.
%   p = 1: uniform in v
%   p = 2: uniform in v^2 (i.e. uniform in energy)
%   p = 3: uniform in v^3 (i.e. uniform in velocity-space volume)

E0       = 1000;                    % eV (bin centre)
Eres_vec = logspace(log10(0.05), 0, 5);     % dE/E0, must be < 2 so that E1 > 0
Eres_vec = logspace(-2, 0, 10);     % dE/E0, must be < 2 so that E1 > 0
ps       = [1 2 3];
N        = 1e6;
v0       = sqrt(2*E0);              % velocity at bin centre (m = 1 units)

% Analytic raw moment E[v^n] when v^p ~ U(a,b), with a = v1^p, b = v2^p
rawmom = @(n,a,b,p) (b.^(n/p+1) - a.^(n/p+1)) ./ ((n/p+1)*(b-a));

nE = numel(Eres_vec);  nP = numel(ps);
vmean_s = zeros(nE,nP);  vvar_s = zeros(nE,nP);   % sampled
vmean_a = zeros(nE,nP);  vvar_a = zeros(nE,nP);   % analytic

for i = 1:nE
  Eres = Eres_vec(i);
  dE = E0*Eres;
  v1 = sqrt(2*(E0 - 0.5*dE));
  v2 = sqrt(2*(E0 + 0.5*dE));

  rng(1);
  r = rand(N,1);                  % same random numbers for all p
  for j = 1:nP
    p = ps(j);
    a = v1^p;  b = v2^p;
    v = (a + r*(b-a)).^(1/p);

    vmean_s(i,j) = mean(v);
    vvar_s(i,j)  = var(v);

    m1 = rawmom(1,a,b,p);
    m2 = rawmom(2,a,b,p);
    vmean_a(i,j) = m1;
    vvar_a(i,j)  = m2 - m1^2;
  end
end

% Use the analytic values for the comparison (no sampling noise),
% and the sampled values as markers to confirm them.
cols = lines(nP);
lbl  = arrayfun(@(p) sprintf('p = %d', p), ps, 'UniformOutput', false);

figure('Position',[100 100 1100 380]);

% 1) Mean velocity relative to the bin-centre velocity
subplot(1,3,1); hold on;
for j = 1:nP
  plot(Eres_vec, vmean_a(:,j)/v0, '-', 'Color', cols(j,:), 'LineWidth', 1.5);
  plot(Eres_vec, vmean_s(:,j)/v0, 'o', 'Color', cols(j,:), 'HandleVisibility','off');
end
set(gca,'XScale','log'); xlabel('E_{res} = \DeltaE/E_0'); ylabel('<v> / v_0');
title('Mean velocity'); legend(lbl,'Location','best'); box on; grid on;

% 2) Spread (sqrt of variance, i.e. the "temperature" proxy) relative to the mean
subplot(1,3,2); hold on;
for j = 1:nP
  plot(Eres_vec, sqrt(vvar_a(:,j))./vmean_a(:,j), '-', 'Color', cols(j,:), 'LineWidth', 1.5);
  plot(Eres_vec, sqrt(vvar_s(:,j))./vmean_s(:,j), 'o', 'Color', cols(j,:), 'HandleVisibility','off');
end
set(gca,'XScale','log','YScale','log'); xlabel('E_{res} = \DeltaE/E_0');
ylabel('\sigma_v / <v>'); title('Relative velocity spread'); box on; grid on;

% 3) Relative difference in variance w.r.t. p = 1
subplot(1,3,3); hold on;
for j = 2:nP
  plot(Eres_vec, abs(vvar_a(:,j)./vvar_a(:,1) - 1)*100, '-', ...
    'Color', cols(j,:), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log'); xlabel('E_{res} = \DeltaE/E_0');
ylabel('|var_p / var_{p=1} - 1|  (%)'); title('Variance error vs. p = 1');
legend(lbl(2:end),'Location','best'); box on; grid on;

% Print the numbers for the resolution in your example
[~,k] = min(abs(Eres_vec - 0.24));
fprintf('E_res = %.2f\n', Eres_vec(k));
for j = 1:nP
  fprintf('p = %d:  <v> = %.3f,  var = %.4f\n', ps(j), vmean_a(k,j), vvar_a(k,j));
end

%%
% Multi-bin test: a Gaussian VDF is binned on a constant-dE/E energy grid,
% then reconstructed by sampling uniformly in v^p inside each bin.
% Compare the recovered mean and spread with the truth as the distribution
% width goes from << bin width to >> bin width.

E0   = 1000;                         % eV, centre of one bin (grid anchor)
Eres = 0.24;                         % dE/E
v0   = sqrt(2*E0);                   % true mean velocity (m = 1 units)
ps   = [1 2 3];

% --- Contiguous grid: bin k spans [Elo(k), Ehi(k)], Ehi(k) = Elo(k+1) ---
r   = (1+Eres/2)/(1-Eres/2);         % ratio of consecutive edges
k   = -40:40;
Elo = E0*(1-Eres/2)*r.^k;
Ehi = Elo*r;                         % bin k=0 is [E0-dE/2, E0+dE/2], as in your code
v1  = sqrt(2*Elo);
v2  = sqrt(2*Ehi);
dv  = v0*Eres/2;                     % bin width in v at v0 (dv/v = dE/2E)

% --- Widths to scan, in units of the bin width ---
sig_rel = logspace(log10(0.05), log10(10), 25);   % sigma_v / dv
vt      = sig_rel*dv;                            % true sigma_v (v0 stays fixed)

Phi     = @(v,mu,s) 0.5*erfc(-(v-mu)/(s*sqrt(2)));        % Gaussian CDF
rawmom  = @(n,a,b,p) (b.^(n/p+1) - a.^(n/p+1)) ./ ((n/p+1)*(b-a));

nS = numel(vt); nP = numel(ps);
mean_a = zeros(nS,nP);  var_a = zeros(nS,nP);

for i = 1:nS
  w = Phi(v2,v0,vt(i)) - Phi(v1,v0,vt(i));     % probability mass per bin
  w = w/sum(w);
  for j = 1:nP
    p  = ps(j);
    a  = v1.^p;  b = v2.^pillusy;
    m1 = rawmom(1,a,b,p);                    % per-bin mean of v
    m2 = rawmom(2,a,b,p);                    % per-bin <v^2>
    mu = sum(w.*m1);
    mean_a(i,j) = mu;
    var_a(i,j)  = sum(w.*m2) - mu^2;         % total variance (mixture)
  end  
  %histogram(a)
  %pause
end

% --- Plots ---
cols = lines(nP);
lbl  = arrayfun(@(p) sprintf('p = %d',p), ps, 'UniformOutput', false);

figure('Position',[100 100 1100 380]);

subplot(1,3,1); hold on;
for j = 1:nP
  plot(sig_rel, (mean_a(:,j)/v0 - 1)*100, 'Color', cols(j,:), 'LineWidth', 1.5);
end
set(gca,'XScale','log'); xlabel('\sigma_v / \Deltav_{bin}');
ylabel('(<v>_{rec} / v_0 - 1)  (%)'); title('Mean velocity error');
legend(lbl,'Location','best'); grid on; box on;

subplot(1,3,2); hold on;
for j = 1:nP
  plot(sig_rel, sqrt(var_a(:,j))./vt(:), 'Color', cols(j,:), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log'); xlabel('\sigma_v / \Deltav_{bin}');
ylabel('\sigma_{rec} / \sigma_{true}'); title('Recovered spread');
grid on; box on;

subplot(1,3,3); hold on;
for j = 2:nP
  plot(sig_rel, (var_a(:,j)./var_a(:,1) - 1)*100, 'Color', cols(j,:), 'LineWidth', 1.5);
end
set(gca,'XScale','log'); xlabel('\sigma_v / \Deltav_{bin}');
ylabel('(var_p / var_{p=1} - 1)  (%)'); title('Effect of p on variance');
legend(lbl(2:end),'Location','best'); grid on; box on;

% --- Sampling check at sigma = 1 bin width (same procedure as your script) ---
[~,i0] = min(abs(sig_rel - 1));
w = Phi(v2,v0,vt(i0)) - Phi(v1,v0,vt(i0));  w = w/sum(w);
cw = [0, cumsum(w)];  cw(end) = 1 + eps;
N = 1e6;  rng(1);
[~,~,idx] = histcounts(rand(N,1), cw);           % pick a bin for each particle
u = rand(N,1);                                   % same uniforms for all p
fprintf('sigma/dv = %.2f\n', sig_rel(i0));
for j = 1:nP
  p = ps(j);
  a = v1(idx).^p;  b = v2(idx).^p;
  v = (a(:) + u.*(b(:)-a(:))).^(1/p);
  fprintf('p=%d: mean %.4f (analytic %.4f),  std %.4f (analytic %.4f)\n', ...
    p, mean(v), mean_a(i0,j), std(v), sqrt(var_a(i0,j)));
end

%%
% Forward-model (bin-integrated Gaussian, Poisson/multinomial ML) vs. moment-based
% estimate of sigma (and hence T ~ sigma^2) on a constant-dE/E energy grid.
% Scan: true width sigma/dv, beam offset relative to the bin grid, finite counts.

rng(1);
E0 = 1000;  Eres = 0.24;
v0 = sqrt(2*E0);

% --- Contiguous grid, bin k=0 is [E0(1-Eres/2), E0(1+Eres/2)] ---
rr  = (1+Eres/2)/(1-Eres/2);
k   = -30:30;
Elo = E0*(1-Eres/2)*rr.^k;
v1  = sqrt(2*Elo);
v2  = sqrt(2*Elo*rr);
i0  = find(k==0);
dv  = v2(i0) - v1(i0);                         % bin width in v at the beam position

% --- Settings ---
sig_rel = logspace(log10(0.03), log10(2), 20); % sigma_true / dv
offsets = [0 0.5];                             % beam centre relative to bin 0 centre, in dv
Ntot    = 2e4;                                 % detected particles per realisation
nrep    = 200;

Phi    = @(v,mu,s) 0.5*erfc(-(v-mu)./(s*sqrt(2)));
binP   = @(mu,s) Phi(v2,mu,s) - Phi(v1,mu,s);  % probability per bin
m1     = (v1+v2)/2;                            % p = 1 sampling inside the bin
m2     = (v1.^2 + v1.*v2 + v2.^2)/3;

% Negative log-likelihood in units of dv; x = [(mu-v0)/dv, log(sigma/dv)]
%nll = @(x,c,pos) -sum(c(pos).*log(max(binP(v0+x(1)*dv, max(exp(x(2)),1e-3)*dv)(pos), realmin)));
opts = optimset('Display','off','TolX',1e-6,'TolFun',1e-9);

nS = numel(sig_rel);  nO = numel(offsets);
R_mom = zeros(nS,nO,3);  R_fit = zeros(nS,nO,3);   % [16 50 84] percentiles of sigma^2 ratio

for io = 1:nO
  for is = 1:nS
    s_true  = sig_rel(is)*dv;
    mu_true = v0 + offsets(io)*dv;
    P       = binP(mu_true, s_true);  P = P/sum(P);
    cw      = [0, cumsum(P)];  cw(end) = 1 + eps;

    rm = zeros(nrep,1);  rf = zeros(nrep,1);
    for n = 1:nrep
      [c,~] = histcounts(rand(Ntot,1), cw);   % multinomial counts
      w   = c/sum(c);  pos = c > 0;

      % (a) moment estimator (reconstructed particles, p = 1)
      mu_m  = sum(w.*m1);
      var_m = sum(w.*m2) - mu_m^2;
      rm(n) = var_m/s_true^2;

      % (b) forward-model ML fit
      x0 = [(mu_m-v0)/dv, log(max(0.5*sqrt(var_m)/dv, 0.05))];
      %xh = fminsearch(@(x) nll(x,c,pos), x0, opts);
      xh = fminsearch(@(x) nll_bins(x, c, pos, binP, v0, dv), x0, opts);
      sh = max(exp(xh(2)),1e-3)*dv;
      rf(n) = sh^2/s_true^2;
    end
    R_mom(is,io,:) = prctile(rm,[16 50 84]);
    R_fit(is,io,:) = prctile(rf,[16 50 84]);
  end
end

% --- Plot ---
figure('Position',[100 100 1000 400]);
cm = [0.84 0.15 0.16];  cf = [0.12 0.47 0.71];
for io = 1:nO
  subplot(1,2,io); hold on;
  fill([sig_rel fliplr(sig_rel)], [R_mom(:,io,1)' fliplr(R_mom(:,io,3)')], cm, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
  fill([sig_rel fliplr(sig_rel)], [R_fit(:,io,1)' fliplr(R_fit(:,io,3)')], cf, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
  plot(sig_rel, R_mom(:,io,2), '-', 'Color', cm, 'LineWidth', 1.5, 'DisplayName','moments');
  plot(sig_rel, R_fit(:,io,2), '-', 'Color', cf, 'LineWidth', 1.5, 'DisplayName','forward fit');
  if offsets(io) == 0       % simple model: sigma_true^2 + dv^2/12
    plot(sig_rel, 1 + 1./(12*sig_rel.^2), 'k--', 'DisplayName','1 + \Deltav^2/(12\sigma^2)');
  end
  yline(1,'k:','HandleVisibility','off');
  set(gca,'XScale','log','YScale','log');
  xlabel('\sigma_{true} / \Deltav_{bin}');
  ylabel('T_{est} / T_{true}  (= \sigma^2 ratio)');
  title(sprintf('beam offset = %.1f bin widths', offsets(io)));
  legend('Location','northeast'); grid on; box on;
end

fprintf('Floor (centred beam): sigma_floor/dv = %.3f -> T_floor/T(sigma=dv) = %.3f\n', ...
  1/sqrt(12), 1/12);


function f = nll_bins(x, c, pos, binP, v0, dv)
mu = v0 + x(1)*dv;
s  = max(exp(x(2)), 1e-3)*dv;
P  = binP(mu, s);                  % temporary variable, no chained ()
Pp = P(pos);
f  = -sum(c(pos) .* log(max(Pp, realmin)));
end