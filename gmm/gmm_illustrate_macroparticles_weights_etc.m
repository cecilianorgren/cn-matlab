totcounts = nansum(iPDist_counts(1).data(:));
MPc = iPDist_counts(1).elim([100 Inf]).macroparticles('ntot',totcounts,'skipzero',1,'counts');
MPf = iPDist(1).elim([100 Inf]).macroparticles('ntot',1e6,'skipzero',1);

hca = subplot(1,1,1);
histogram(hca,MPf.vx,-3000:200:3000,'Normalization','pdf')
hold(hca,'on')
histogram(hca,MPc.vx,-3000:200:3000,'Normalization','pdf')
hold(hca,'off')
legend({'Equal dn for all MPs','1 count per MP'},'box','off')
hca.XLabel.String = 'v_x (km/s)';
hca.YLabel.String = 'Probability distribution function';

%%
nMP = [1e5 2e5 5e5 1e6 2e6 5e6];
clear MPf;
for iMP = 1:numel(nMP)
  MPf{iMP} = iPDist(1).elim([100 Inf]).macroparticles('ntot',nMP(iMP),'skipzero',1);
end

dn = cellfun(@(x) x.df.*x.dv, MPf, 'UniformOutput', false);

NMP = numel(nMP);
h = setup_subplots(NMP,1);
isub = 1;
for iMP = 1:NMP
  hca = h(isub); isub = isub + 1;
  histogram(hca,dn{iMP},20,'Normalization','cdf')
  hca.XLabel.String = 'dn=dfv';
  hca.YLabel.String = 'CDF';
  irf_legend(hca,sprintf('n_{MP} = %g',nMP(iMP)),[0.02 0.7],'k')
end

