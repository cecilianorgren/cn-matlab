% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

%%
events = struct([]);
iEvent = 0;
if 1 % Torbert 2018
  iEvent = iEvent + 1;  
  tint = irf.tint('2017-07-11T22:33:45.00Z/2017-07-11T22:34:20.00Z'); 
  L = [0.9482,-0.255,-0.1893];
  M = [0.1818,0.9245,-0.3350];
  N = [0.2604,0.2832,0.9239];
  lmn = [L; M; N];

  ic = 1:4;
  ve = cell(1,numel(ic));
  for mms_id = 1:numel(ic)
    ve{mms_id} = mms.get_data('Ve_gse_fpi_brst_l2',tint,mms_id);
    Pe{mms_id} = mms.get_data('Pe_gse_fpi_brst_l2',tint,mms_id);
    Te{mms_id} = mms.get_data('Te_gse_fpi_brst_l2',tint,mms_id);
    B{mms_id}  = mms.get_data('B_gse_brst_l2',tint,mms_id);
    facPepp{mms_id} = mms.rotate_tensor(Pe{mms_id},'fac',B{mms_id},'pp'); % Peperp1 = Peperp2
    facTepp{mms_id} = mms.rotate_tensor(Te{mms_id},'fac',B{mms_id},'pp'); % Peperp1 = Peperp2
    Q_{mms_id} = (facPepp{mms_id}.xy.data.^2+facPepp{mms_id}.xz.data.^2+facPepp{mms_id}.yz.data.^2)./(facPepp{mms_id}.yy.data.^2+2*facPepp{mms_id}.yy.data.*facPepp{mms_id}.xx.data);
    Q{mms_id} = irf.ts_scalar(ve{mms_id}.time,Q_{mms_id});
  end

  events(iEvent).tint = tint;
  events(iEvent).mms_id = ic;
  events(iEvent).lmn = lmn;
  events(iEvent).ve = ve;
  events(iEvent).facPepp = facPepp;
  events(iEvent).facTepp = facTepp;
  events(iEvent).sqrtQ = Q;  
end
if 0 % Ergun 2022
  iEvent = iEvent + 1;  
  tint = irf.tint('2018-08-27T12:15:35.00Z/2018-08-27T12:15:55.00Z'); 
  L = [0.910,-0.385,-0.155];
  M = [0.415,0.848,0.331];
  N = [0.004,-0.365,0.931];
  lmn = [L; M; N];

  ic = 1:4;
  ve = cell(1,numel(ic));
  for mms_id = 1:numel(ic)
    ve{mms_id} = mms.get_dat1a('Ve_gse_fpi_brst_l2',tint,mms_id);
    Pe{mms_id} = mms.get_data('Pe_gse_fpi_brst_l2',tint,mms_id);
    B{mms_id}  = mms.get_data('B_gse_brst_l2',tint,mms_id);
    facPepp{mms_id} = mms.rotate_tensor(Pe{mms_id},'fac',B{mms_id},'pp'); % Peperp1 = Peperp2
    Q_{mms_id} = (facPepp{mms_id}.xy.data.^2+facPepp{mms_id}.xz.data.^2+facPepp{mms_id}.yz.data.^2)./(facPepp{mms_id}.yy.data.^2+2*facPepp{mms_id}.yy.data.*facPepp{mms_id}.xx.data);
    Q{mms_id} = irf.ts_scalar(ve{mms_id}.time,Q_{mms_id});
  end

  events(iEvent).tint = tint;
  events(iEvent).mms_id = ic;
  events(iEvent).lmn = lmn;
  events(iEvent).ve = ve;
  events(iEvent).sqrtQ = Q;  
end

% Add others...


%%
units = irf_units;
fontsize = 16;

nEvents = numel(events);
h = setup_subplots(2,1);
isub = 1;

if 0 % Scatter Q(vL,vM)
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0,'+','markersize',20,'linewidth',3,'color','k')
  hold(hca,'on')
  for iEvent = 1:nEvents
    event = events(iEvent);
    for ic = 1:numel(event.mms_id)
      nSmooth = 5;
      ve = event.ve{ic}.smooth(nSmooth)*event.lmn';
      Q = event.sqrtQ{ic}.smooth(nSmooth);
      scatter(hca,ve.x.data,ve.y.data,20,Q.data,'filled')
    end
  end
  
  hcb = colorbar(hca);
  hcb.YLabel.String = '$\sqrt{Q}$';
  hcb.YLabel.Interpreter = 'latex';
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'off')
  axis(hca,'equal')
  axis(hca,'square')
  colormap(pic_colors('candy_gray'))
  hca.XTick = -20e3:5e3:20e3;
  hca.YTick = -20e3:5e3:20e3;
  hca.XLabel.String = 'v_{eL} (km/s)';
  hca.YLabel.String = 'v_{eM} (km/s)';
end
if 1 % Scatter vLM(vL,vM)
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0,'+','markersize',20,'linewidth',3,'color','k')
  hold(hca,'on')
  for iEvent = 1:nEvents
    event = events(iEvent);
    for ic = 1:numel(event.mms_id)
      nSmooth = 5;
      ve = event.ve{ic}.smooth(nSmooth)*event.lmn';
      ve_LM = sqrt(ve.x.data.^2 + ve.y.data.^2);
      scatter(hca,ve.x.data,ve.y.data,20,ve_LM,'filled')
    end
  end
  
  hcb = colorbar(hca);
  hcb.YLabel.String = '$\sqrt{v_{eL}^2 + v_{eM}^2}$';
  hcb.YLabel.Interpreter = 'latex';
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'off')
  axis(hca,'equal')
  axis(hca,'square')
  colormap(pic_colors('candy_gray'))
  hca.XTick = -20e3:5e3:20e3;
  hca.YTick = -20e3:5e3:20e3;
  hca.XLabel.String = 'v_{eL} (km/s)';
  hca.YLabel.String = 'v_{eM} (km/s)';
end
if 1 % Scatter Teperp(vL,vM)
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0,'+','markersize',20,'linewidth',3,'color','k')
  hold(hca,'on')
  for iEvent = 1:nEvents
    event = events(iEvent);
    for ic = 1:numel(event.mms_id)
      nSmooth = 5;
      Te = event.facTepp{ic}.yy.smooth(nSmooth).data; % Teperp1 = Teperp2
      vt = sqrt(2*units.eV*Te/units.me)*1e-3;
      scatter(hca,ve.x.data,ve.y.data,20,vt,'filled')
    end
  end
  
  hcb = colorbar(hca);
  hcb.YLabel.String = 'v_{te\perp}';
  hcb.YLabel.Interpreter = 'tex';
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'off')
  axis(hca,'equal')
  axis(hca,'square')
  colormap(pic_colors('candy_gray'))
  hca.XTick = -20e3:5e3:20e3;
  hca.YTick = -20e3:5e3:20e3;
  hca.XLabel.String = 'v_{eL} (km/s)';
  hca.YLabel.String = 'v_{eM} (km/s)';
end
if 0 % angle
  hca = h(isub); isub = isub + 1;
  plot(hca,0,0,'+','markersize',20,'linewidth',3,'color','k')
  hold(hca,'on')
  for iEvent = 1:nEvents
    event = events(iEvent);
    for ic = 1:numel(event.mms_id)
      nSmooth = 1;
      ve = event.ve{ic}.smooth(nSmooth)*event.lmn';
      ve_angle = acosd(ve.x.data./sqrt(ve.x.data.^2 + ve.y.data.^2));
      Q = event.sqrtQ{ic}.smooth(nSmooth);
      %scatter(hca,ve.x.data,ve.y.data,5,ve_angle)
      plot(hca,ve_angle,Q.data,'*')
    end
  end
  
  %hcb = colorbar(hca);
  %hcb.YLabel.String = '\sqrt{Q}';
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'off')
  axis(hca,'equal')
  axis(hca,'square')
  colormap(pic_colors('candy_gray'))
  
end


c_eval('h(?).FontSize = fontsize;',1:numel(h))