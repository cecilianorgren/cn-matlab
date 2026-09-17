%[Drms,partitions,part_str] = gmm_gaussianity_loss_all_partitions(gm(65:66,2));

Nparts = numel(partitions{1});
[~,isort] = sort(Drms{1});
nGroups = cellfun(@(x)numel(x),partitions{1});
[a,p1,p2] = plotyy(1:Nparts,Drms{1}(isort),1:Nparts,nGroups(isort));
p1.Marker = '*';
p2.Marker = '*';
a(1).XLim = [0 Nparts+1];
a(2).XLim = [0 Nparts+1];
a(1).YLim = [0 0.6];
a(2).YLim = [0 6];
a(1).FontSize = 15;
a(2).FontSize = 15;
xticks(1:Nparts);
xticklabels(part_str{1}(isort))
grid on
yticks(a(1),0:0.1:0.6)
yticks(a(2),0:1:6)
a(1).YLabel.String = 'D RMS';
a(2).YLabel.String = 'Number of groups';
a(1).XLabel.String = 'Partition of components into groups';