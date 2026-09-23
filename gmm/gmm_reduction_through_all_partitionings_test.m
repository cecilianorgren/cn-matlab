tic
iK = 5;
[Drms,partitions,part_str] = gmm_gaussianity_loss_all_partitions(gm(:,iK));
toc

nt = numel(Drms);
nPartitions = numel(Drms{1});
nGroups_all = zeros(nt,nPartitions);
D_all = zeros(nt,nPartitions);

for it = 1:numel(Drms)
  Nparts = numel(partitions{1});
  [~,isort] = sort(Drms{it});  
  nGroups = cellfun(@(x)numel(x),partitions{1});
  nGroups_all(it,:) = nGroups(isort);
  D_all(it,:) = Drms{it}(isort);
  if 1 % plot
    [a,p1,p2] = plotyy(1:Nparts,Drms{it}(isort),1:Nparts,nGroups(isort));
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
    %pause
  end
end
%%
h = irf_plot(6);


if 1 % B
  hca = irf_panel('B');
  irf_plot(hca,gseB)
  set(hca,'ColorOrder',mms_colors('xyza'))
  hca.YLabel.String = 'B (nT)';
end
if 1 % dEFlux ion PD2
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
end
if 1 % BIC(K)
  hca = irf_panel('BIC');
  specrec.p = cat(2,DB.events(id_event).runs(:).BIC);
  specrec.f = cat(2,DB.events(id_event).runs(:).NumComponents);  
  specrec.t = PD.time.epochUnix;
  irf_spectrogram(hca,specrec,'lin')  
  hca.YLabel.String = 'K';
  hcb = colorbar(hca);
  hcb.YLabel.String = 'BIC';
end
if 1 % rel L2 error 
  hca = irf_panel('Rel L2 error');
  specrec.p = D_all;
  specrec.f = 1:numel(partitions{1});  
  specrec.t = PD.time.epochUnix;
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel')  
  hca.YLabel.String = 'Partition id';
  hcb = colorbar(hca);
  hcb.YLabel.String = 'Rel L2 error';
end
if 1 % rel L2 error , corresponding n of groups
  hca = irf_panel('Rel L2 error, ngroups');
  specrec.p = nGroups_all;
  specrec.f = 1:numel(partitions{1});  
  specrec.t = PD.time.epochUnix;
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel')  
  hca.YLabel.String = 'Partition id';
  hcb = colorbar(hca);
  hcb.YLabel.String = 'N groups';
end
if 1 % rel L2 error , corresponding n of groups, masked
  hca = irf_panel('Rel L2 error, ngroups, masked');
  threshold = 0.2;
  specrec.p = nGroups_all;
  specrec.p(D_all>threshold) = NaN;
  specrec.f = 1:numel(partitions{1});  
  specrec.t = PD.time.epochUnix;
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel')  
  hca.YLabel.String = 'Partition id';
  hcb = colorbar(hca);
  hcb.YLabel.String = 'N groups';
  irf_legend(hca,sprintf('Rel L2 error<%g',threshold),[0.98 0.98],'color','k')
end

c_eval('h(?).Layer = ''top'';',1:numel(h))
colormap([flipdim(irf_colormap('Spectral'),1)])
irf_plot_axis_align(h)
irf_zoom(h,'x',PD.time([1 PD.length]))
h(end).XTickLabelRotation = 0;