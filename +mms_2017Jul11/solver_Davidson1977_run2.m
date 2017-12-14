% run for array of Te and E, 3rd dim is k
n = 0.06e6;
ne = n;
ni = n;
Te_vec = 500:100:700; Te_K_vec = Te_vec*units.eV/units.kB;  % use perpendicular temperature -> diamagnetic drift
Ti_vec = 7000; Ti_K_vec = Ti_vec*units.eV/units.kB; % use perpendicular temperature -> diamagnetic drift
B = 10e-9; 

E_vec = (-50:10:-10)*1e-3;
LnLT = 0.1; % Ln/LT, Ln generally smaller then LT

nTe = numel(Te_vec);
nE = numel(E_vec);

nk = 200;

wr_allmax = nan(nE,nTe);
wi_allmax = nan(nE,nTe);
k_allmax = nan(nE,nTe);
vph_allmax = nan(nE,nTe);

wr_all = nan(nE,nTe,nk);
wi_all = nan(nE,nTe,nk);
k_all = nan(nE,nTe,nk);
vph_all = nan(nE,nTe,nk);

for iE = 1:nE
  for iTe = 1:nTe
    Te = Te_vec(iTe);
    Te_K = Te_K_vec(iTe);
    E = E_vec(iE);
    mms_2017Jul11.solver_Davidson1977
    wr_allmax(iE,iTe) = wrmax;
    wi_allmax(iE,iTe) = wimax;
    k_allmax(iE,iTe) = kmax;
    vph_allmax(iE,iTe) = vphmax;
    
    wr_all(iE,iTe,:) = wr_store;
    wi_all(iE,iTe,:) = wi_store;
    k_all(iE,iTe,:) = kmax;
    vph_all(iE,iTe,:) = vphmax;
  end
end

roe_all = sqrt(2*qe*Te_vec/me)/(qe*B/me);
roe_mat = repmat(roe_all,nE,1);