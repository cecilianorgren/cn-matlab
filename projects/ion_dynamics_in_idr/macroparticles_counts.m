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

%%
it = 1000;
pdist = iPDist(it);
counts = iPDist_counts(it);
totCounts = nansum(counts.data(:));
MP = pdist.macroparticles;
MPc = counts.macroparticles;

%unique(counts.data(not(isnan(counts.data))))

%counts.data(isnan(counts.data)) = 0; % needs to be implemented inside PDist, done
MPc_counts_random = counts.macroparticles('ntot',totCounts,'ntot_division','counts','positioning','random','skipzero',1);
MPc_counts_center = counts.macroparticles('ntot',totCounts,'ntot_division','counts','positioning','center','skipzero',1);
unique(MPc_counts_random.df) % I thought this would give 1, so something to double check
unique(MPc_counts_center.df) % I thought this would give 1, so something to double check

%%


%MPc_counts_random = counts.macroparticles('ntot',totCounts,'ntot_division','counts','positioning','random','skipzero',1);

% for counts, MP_c_c.df are expected to be integers (or close enough to, as given by the counts input)
MP_c_c = counts.macroparticles('ntot',totCounts,'ntot_division','counts','positioning','center','skipzero',1)
MP_c_r = counts.macroparticles('ntot',totCounts,'ntot_division','counts','positioning','random','skipzero',1)

MP_f_c = pdist.macroparticles('ntot',totCounts,'ntot_division','fdv','positioning','center','skipzero',1)
MP_f_r = pdist.macroparticles('ntot',totCounts,'ntot_division','fdv','positioning','random','skipzero',1)

% There is some discrepancy here, typically it was often a factor 2, now
% larger, maybe PDist.d3v (used to get dv, phase space volume of each particle) is wrong?
n_c = sum(MP_f_c.df.*MP_f_r.dv)*1e-6
n_r = sum(MP_f_r.df.*MP_f_r.dv)*1e-6

ne(it).data

hca = subplot(1,3,1);
plot(hca,MP_c_c.vx,MP_c_r.vx,'.')
hold(hca,'on')
plot(hca,hca.XLim,hca.XLim)
hold(hca,'off')

hca = subplot(1,3,2);
scatter(hca,MP_c_c.vx,MP_c_c.vy,'.')

hca = subplot(1,3,3);
scatter(hca,MP_c_r.vx,MP_c_r.vy,'.')

h = findobj(gcf,'type','axes'); h = h(end:-1:1);

linkprop(h(2:3),{'XLim','YLim'})
