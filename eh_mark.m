% make multiple markings

figure; h=irf_plot({diE3(1:fix(end/3),:),diE4(1:fix(end/3),:)},'comp');

[tint_eh nscobs comments] = eh_tint; % load time intervals
ax=h; % plot to be marked is
% seen on one/seen on 2/developing seen on one/developing seen on two
color={'green',[1 1 0],[0 0.4 1],'red','k','k','k'}; 
for k = 1:size(tint_eh,1)
    color{nscobs(k)};
    irf_pl_mark(ax,tint_eh{k},color{nscobs(k)})
end