% 
PD = iPDist3;
PD_counts = iPDist3_counts;
elim_lowE = [0 200];
av_counts = mean(sum(PD_counts.nan2zero.data,[3 4]),1);
counts_25 = prctile(sum(PD_counts.nan2zero.data,[3 4]),25,1);
counts_50 = prctile(sum(PD_counts.nan2zero.data,[3 4]),50,1);
counts_75 = prctile(sum(PD_counts.nan2zero.data,[3 4]),75,1);

counts_E = sum(PD_counts.nan2zero.data(:,:,:,:),[3 4]);

av_counts_lowE = mean(sum(PD_counts.elim(elim_lowE).nan2zero.data,[3 4]),1);
av_counts_lowE_1 = mean(av_counts_lowE);

N_rem_per_elevel = round(av_counts_lowE_1);

f_per_count = PD.data./PD_counts.data;

f_per_count_1 = nanmean(nanmean(f_per_count,[3 4]),1);

%%
PD_new = PD;
for it = 1:PD.length
  it
  for iE = 1:32
    if sum(PD_counts.data(it,iE,:)) == N_rem_per_elevel

    nRem = N_rem_per_elevel;
    iRem = 0;
    while iRem < nRem
      iAz = randi(32,1);
      iPol = randi(16,1);
      if PD.data(it,iE,iAz,iPol) == 0
        continue
      else
        PD_new.data(it,iE,iAz,iPol) = PD_new.data(it,iE,iAz,iPol) - f_per_count(it,iE,iAz,iPol);
        iRem = iRem + 1;
      end
    end
  end
end

%%
fontsize = 14;

hca = subplot(2,1,1);
plot(hca,PD_counts.depend{1}(1,:),av_counts)
hold(hca,'on')
hl = plot(hca,PD_counts.elim(elim_lowE).depend{1}(1,:),av_counts_lowE,'*');
hold(hca,'off')
hca.XLabel.String = 'Energy (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = 'Counts summed over angles and averaged over time';
hca.XScale = 'log';
irf_legend(hca,{sprintf('Average at low energies: %.2f',av_counts_lowE_1)},[0.02 0.98],'color',hl.Color,'fontsize',fontsize)
hca.FontSize = fontsize;


hca = subplot(2,1,2);
plot(hca,PD_counts.depend{1}(1,:),[counts_25;counts_50;counts_75])
%hold(hca,'on')
%hl = plot(hca,PD_counts.elim(elim_lowE).depend{1}(1,:),av_counts_lowE,'*');
%hold(hca,'off')
hca.XLabel.String = 'Energy (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = 'Counts summed over angles and percentiled over time';
hca.XScale = 'log';
%irf_legend(hca,{sprintf('Average at low energies: %.2f',av_counts_lowE_1)},[0.02 0.98],'color',hl.Color,'fontsize',fontsize)
hca.FontSize = fontsize;