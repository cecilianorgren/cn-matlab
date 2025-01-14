%% Load data
mms.db_init('local_file_db','/Volumes/mms');

tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); 

ic = 3;
ePDist = mms.get_data('PDe_fpi_brst_l2',tint,ic);
ne = mms.get_data('Ne_fpi_brst_l2',tint,ic);

iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
iPDist_counts = iPDist; iPDist_counts.data = (iPDist_counts.data./iPDistErr.data).^2;

c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

%% Make macroparticles
% Number of macroparticles, although inside PDist, the minimum is set to 
% one per bin woth non-zero f. The more you have, the more evenly they will 
% be weighted, see hist(mp.dn) below.
nMP = 1e5; 

% Pick your distribution(s), PDist.macroparticles returns an 1 x nt struct array
it = 1200; % it = 1200:1201;
pdist = ePDist(it);

% Generate macroparticles
% It might be necessary to remove the lowest energy levels due to
% instrumental contamination
elim = [100 Inf];
pdist = pdist.elim(elim);
mp = pdist.macroparticles('ntot',nMP,'scpot',scPot(it));
% mp.DepX shows which FPI bin they were generated from.

% Compare sum of partial densities with FPI density
disp(sprintf('ne_fpi = %g, sum(dn_mp) = %g',ne(it).data,sum(mp.dn)))

% Now, do what you want with the macroparticles
VX = mp.vx;
VY = mp.vy;
VZ = mp.vz;

% Rotate into a new coordinate system, for example a field aligned system
L = [0 1 0];
M = [0 0 1];
N = [1 0 0];
R = [L; M ; N];

% Rotate
VX2 = [VX VY VZ]*L';
VY2 = [VX VY VZ]*M';
VZ2 = [VX VY VZ]*N';

% Same thing
VXYZ = [VX VY VZ]*R';

VABS = sqrt(VX.^2 + VY.^2 + VZ.^2);

% Bin macroparticles into new grid
% Same grid as FPI, here just using the bin numbers directly
nDep1 = numel(pdist.depend{1}(1,:));
nDep2 = numel(pdist.depend{2}(1,:));
nDep3 = numel(pdist.depend{3}(1,:));
[sum_df_fpi,~,mid,~] = histcn([mp.iDep1, mp.iDep2, mp.iDep3], 0.5:1:(nDep1+0.5), 0.5:1:(nDep2+0.5), 0.5:1:(nDep3+0.5), 'AccumData', mp.df, 'Fun', @sum);
% sum_df_fpi should just give back the original f for each bin: compare
iAz = 12;
iPol = 12;
[squeeze(pdist.data(1,:,iAz,iPol))' sum_df_fpi(:,iAz,iPol)]' % compare f(E) for one of the angular bins

% Cartesian grid
dv = 2000; % km/s
vmax = 100000;
v_edges = -vmax:dv:vmax;
% density contribution (partial density) for each new bin
[sum_dn,~,v_mid,~] = histcn([VX VY VZ], v_edges, v_edges, v_edges, 'AccumData', mp.dn, 'Fun', @sum); % see >> help histcn

% There is something wrong with the units (see comparison of 2D reduced distributions below), factor 10^6?
d3v_kms = dv*3;
d3v_cms = d3v_kms*(1e5)^3; % km^3 -> cm^3
f_cart = sum_dn/d3v_cms; % #/(cm^6/s^3) 
dv_cms = dv*1e5;
f_cart_2d = squeeze(sum(f_cart,3))*dv_cms; 

%% Figure
nRows = 3;
nCols = 2;
isub = 1;

step = 100;

hca = subplot(nRows,nCols,isub); isub = isub + 1;
scatter3(hca,VX(1:step:end),VY(1:step:end),VZ(1:step:end),5,mp.iDep1(1:step:end)) % energy index (Dep1)
hca.XLim = prctile(VABS(:),90)*[-1 1];
hca.YLim = prctile(VABS(:),90)*[-1 1];
hca.ZLim = prctile(VABS(:),90)*[-1 1];

hca = subplot(nRows,nCols,isub); isub = isub + 1;
scatter3(hca,VX(1:step:end),VY(1:step:end),VZ(1:step:end),5,mp.iDep2(1:step:end)) % azimuthal angle index (Dep2)
hca.XLim = prctile(VABS(:),90)*[-1 1];
hca.YLim = prctile(VABS(:),90)*[-1 1];
hca.ZLim = prctile(VABS(:),90)*[-1 1];

hca = subplot(nRows,nCols,isub); isub = isub + 1;
scatter3(hca,VX(1:step:end),VY(1:step:end),VZ(1:step:end),5,mp.iDep3(1:step:end)) % polar angle index (Dep3)
hca.XLim = prctile(VABS(:),90)*[-1 1];
hca.YLim = prctile(VABS(:),90)*[-1 1];
hca.ZLim = prctile(VABS(:),90)*[-1 1];

hca = subplot(nRows,nCols,isub); isub = isub + 1;
hist(hca,mp.dn,100)
hca.XLabel.String = 'dn (= macroparticle weight)';
hca.YLabel.String = 'Number of macro particles';

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot_data = f_cart_2d';
plot_data(plot_data==0) = NaN;
plot_data = plot_data*10^10; % 1/cm^5 -> 1/m^5
pcolor(hca,v_mid{1}*1e-3,v_mid{2}*1e-3,log10(plot_data))
shading(hca,'flat')
hcb = colorbar(hca);
hca.XLabel.String = 'v (10^3 km/s)';
hca.YLabel.String = 'v (10^3 km/s)';
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.Title.String = 'PDist.macroparticles + cart. binning';

hca = subplot(nRows,nCols,isub); isub = isub + 1;
vdf_2d = pdist.reduce('2D',[1 0 0],[0 1 0]);
vdf_2d.plot_plane(hca)
hca.Title.String = 'PDist.reduce';
hca.XLim = 0.99*vmax*1e-3*[-1 1];
hca.YLim = 0.99*vmax*1e-3*[-1 1];


% hca = subplot(nRows,nCols,isub); isub = isub + 1;
% hca = subplot(nRows,nCols,isub); isub = isub + 1;

