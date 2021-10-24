function out = df_selection(df_data,criteria)
% df_selection
% make some kind of selection
t = [df_data.it];
all_opt_best = [df_data.opt_best];
opt_best_B0 = all_opt_best(1,:); 
opt_best_B1 = all_opt_best(2,:); 
opt_best_dt = all_opt_best(3,:); 
opt_xcorr = [df_data.opt_xcorr];

selection = zeros(size(t));
% ind = opt_best_dt > criteria{1}(1) & ... % dt
%       opt_best_dt < criteria{1}(2) & ...
%       opt_best_B0 > criteria{2}(1) & ... % B0
%       opt_best_B0 < criteria{2}(2) & ...
%       opt_best_B1 > criteria{3}(1) & ... % B1
%       opt_best_B1 < criteria{3}(2) & ...
%       opt_xcorr   > criteria{4}(1) & ... % xcorr
%       opt_xcorr   < criteria{4}(2);

ind = opt_best_dt*2           > criteria{1}(1) & ... % 2*dt
      opt_best_dt*2           < criteria{1}(2) & ...
      opt_best_B0-opt_best_B1 > criteria{2}(1) & ... % B0-B1, lowest B
      opt_best_B0-opt_best_B1 < criteria{2}(2) & ...
      opt_best_B1*2           > criteria{3}(1) & ... % 2*B1, change in B across front
      opt_best_B1*2           < criteria{3}(2) & ...
      opt_xcorr               > criteria{4}(1) & ... % xcorr
      opt_xcorr               < criteria{4}(2);
selection(ind) = 1;

% find consecutive sequences and find larges correlation within each of
% them
[B, N, BI] = RunLength(selection);
% find all with slection==1
idf = [];
all_seq_sel = find(B==1);
for isec = 1:numel(B)
  if B(isec) == 1
    part_data = opt_xcorr(BI(isec):(BI(isec)+N(isec)-1));
    [maxval,i_max_corr] = max(part_data);
    idf(end+1) = BI(isec)+i_max_corr-1;
  else
    continue, 
  end
end

if 0 % plot
  plind = 1:numel(opt_xc0rr);
  nrows = 4;
  hca = subplot(nrows,1,1)
  plot(hca,plind,opt_xcorr)
  
  hca = subplot(nrows,1,1)
  plot(hca,plind,selection)
  hold(hca,'on')
  for idf = 1:numel(idf)
    plot(hca,idf,1)
  end
  hold(hca,'off')
end
% Save all data in some readable way
df_save = df_data(idf);

out = df_save;
