energy = ePDist1.ancillary.energy(1,:);
delta_energy_minus = ePDist1.ancillary.delta_energy_minus(1,:);
delta_energy_plus = ePDist1.ancillary.delta_energy_plus(1,:);

energy_minus = energy - delta_energy_minus;
energy_plus = energy + delta_energy_plus;
delta_energy = (delta_energy_plus+delta_energy_minus);

energy_bins = 1:numel(energy);

h = setup_subplots(2,2);
isub = 1;

if 0
hca = h(isub); isub = isub + 1;
[ax,h1,h2] = plotyy(hca,energy_bins,energy*1e-3,energy_bins,(energy))
hca.Title.String = 'Energy (E) of energy channels';
hca.XLabel.String = 'Energy bin';
ax(1).YLabel.String = 'E (keV)';
ax(2).YLabel.String = 'E (eV)';
ax(2).YScale = 'log';
ax(2).YTick = 10.^[0:4];
end

hca = h(isub); isub = isub + 1;
[ax,h1,h2] = plotyy(hca,energy_bins,[energy_minus; energy; energy_plus]*1e-3,energy_bins,[energy_minus; energy; energy_plus]);
hca.Title.String = 'Energy ranges of energy channels';
hca.XLabel.String = 'Energy bin';
hca.XLim = [0 numel(energy)];
ax(1).YLabel.String = '\Delta E (keV)';
ax(2).YLabel.String = '\Delta E (eV)';
ax(2).YScale = 'log';
ax(2).YTick = 10.^[0:4];
legend(hca,{'E-E^-','E','E+E^+','E-E^-','E','E+E^+'},'location','northwest','box','off')

hca = h(isub); isub = isub + 1;
hl = plot(hca,energy_minus(2:end)*1e-3,energy_plus(1:end-1)*1e-3);
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'E^-(2:32) (keV)';
hca.YLabel.String = 'E^+(1:32) (keV)';
%hca.YLim = [0 32];
hca.Title.String = 'Energy ranges of energy channels';
irf_legend(hca,{'The upper energy of the lower channel (y-axis)';'is the same as the lower energy';'of the upper channel (x-axis)'},[0.02 0.98],'color','k')

hca = h(isub); isub = isub + 1;
[ax,h1,h2] = plotyy(hca,energy_bins,delta_energy*1e-3,energy_bins,(delta_energy));
hca.Title.String = 'Energy width (\Delta E=E^++E^-) of channels';
hca.XLabel.String = 'Energy bin';
ax(1).YLabel.String = '\Delta E (keV)';
ax(2).YLabel.String = '\Delta E (eV)';
ax(2).YScale = 'log';
ax(2).YTick = 10.^[0:4];
hca.XLim = [0 numel(energy)];

hca = h(isub); isub = isub + 1;
plot(hca,energy_bins,delta_energy./energy)
hca.Title.String = '\Delta E/E';
hca.YLabel.String = '\Delta E/E';
hca.XLabel.String = 'Energy bin';
hca.XLim = [0 numel(energy)];
irf_legend(hca,{'Except for some oscillations at lower energy channels,';'the energy width over the central energy is constant'},[0.02 0.98],'color','k')

for ip = 1:numel(h);
  h(ip).FontSize = 12;
  %h(ip).XLim = [0 numel(energy)];
end