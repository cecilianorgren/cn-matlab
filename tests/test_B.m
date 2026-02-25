%% load data
[B3offs B4offs]=av_combine('1','offs');
[B3nooffs B4nooffs]=av_combine('1','nooffs');
%% comp
figure;h=irf_plot({B3offs,B4offs,B3nooffs,B4nooffs},'comp')
legend(h(1),'offs C3','offs C4','no offs C3','no offs C4')
title(h(1),'B GSM')
%% comp fgm
figure;irf_plot({gsmB3fgm,gsmB4fgm},'comp');legend('C3','C4')
%%
dBa=irf_add(1,B3,-1,B4);
dBb=irf_abs(irf_add(1,gsmB3,-1,gsmB4));
%%
figure;h=irf_plot({dBa,dBb},'comp')
legend(h(2),'script','c fgm staff combine')
%%
figure;h=irf_plot({B3,B4,gsmB3,gsmB4},'comp')
legend(h(1),'script C3','script C4','c fgm staff combine C3','c fgm staff combine C4')
title(h(1),'B GSM')