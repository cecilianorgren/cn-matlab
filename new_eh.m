%% Define time interval to study
tint=[toepoch(t1) toepoch(t2)];

%% Compare two satellites
% Assumes tint, facE3, facE4 exists
set(gca,'ColorOrder',[[0 0.7 0];[0 0 1]]);
irf_plot({irf_tlim(facE3(:,[1 4]),tint(1),tint(2)),irf_tlim(facE4(:,[1 4]),tint(1),tint(2))},'comp')
irf_zoom(gca,'x',tint);
grid(gca,'off')
set(gca,'fontsize',14);
irf_legend(gca,{'C3','C4'},[0.02 0.1],'fontsize',16)
ylabel('E_{||} [mV/m]','fontsize',14);
xlabel('31-Aug-2007','fontsize',14)
set(gcf,'PaperPositionMode','auto');
%%
print -dpng /Users/Cecilia/EH/EPS-ICPP/epar_comp_3.png