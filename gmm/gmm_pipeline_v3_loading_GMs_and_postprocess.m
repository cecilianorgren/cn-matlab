%% Load DB table
localuser = 'cecilianorgren';
%file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2022.nc'];


%% Loop events and plot brief overview
units = irf_units;

mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Volumes/mms');

doPrint = 0;
doPlot = 1;

% Turn off irf._log warning messages
irf.log('off')

% Save path
pathSave = ['/Users/' localuser '/IRFU Dropbox/Cecilia Norgren/Data/GMM_v1/'];
listFiles = dir([pathSave '*iFF*']);
nFiles = numel(listFiles);

for iFile = 1%:nFiles
  fileName = [listFiles(iFile).folder filesep listFiles(iFile).name];
  GM = load(fileName,'GM_save').GM_save; 

  GBICmin = GM.selectCol(GM.kneeK('method','min'));
  GBICchord = GM.selectCol(GM.kneeK('method','chord'));  
  GMaxK = GM.selectCol(max(GM.nComp,[],'all'));
  
  GBICminCold = GBICchord.subMixture(GBICchord.cold('maxColdEV',1000));
  GBICchordCold = GBICchord.subMixture(GBICchord.cold('maxColdEV',1000));
  GMaxKCold = GMaxK.subMixture(GMaxK.cold('maxColdEV',1000));

  Gm = GBICchord.subMixture(GBICchord.cold('maxColdEV',1000));

  if doPlot
    %% Plot, BIC
    figure(201)
    str_print = [listFiles(iFile).name, '_BIC_overview'];
    h = irf_plot(6);
  
    hca = irf_panel('omni fpi');
    irf_spectrogram(hca,GM.userdata.pdist_derived.omni.deflux.specrec)
  
    hca = irf_panel('fz fpi');
    irf_spectrogram(hca,GM.userdata.pdist_derived.fz.specrec)
  
    if 0
    hca = irf_panel('BIC/maxBIC)');
    irf_spectrogram(hca,GM.tsBIC('normalization','max').specrec( ...
    'p_label','log BIC/max(BIC)','f_label','K'),'log')
    hold(hca,'on')
    irf_plot(hca,GM.tsBICKnee,'k')
    hold(hca,'off')
    end
    if 1
    hca = irf_panel('BIC/minBIC');
    irf_spectrogram(hca,GM.tsBIC('normalization','min').specrec( ...
    'p_label','log BIC-min(BIC)','f_label','K'),'log')
    hold(hca,'on')
    irf_plot(hca,GM.tsBICKnee('method','min'),'color','k','linewidth',0.5)
    irf_plot(hca,GM.tsBICKnee('method','chord'),'color','r','linewidth',0.5)
    hold(hca,'off')
    end
  
    hca = irf_panel('fz by BIC min');
    irf_spectrogram(hca,GBICmin.tsReduced(3).specrec,'log')
  
    hca = irf_panel('fz by BIC chord');
    irf_spectrogram(hca,GBICchord.tsReduced(3).specrec,'log')
  
    hca = irf_panel('fz by max K');  
    irf_spectrogram(hca,GMaxK.tsReduced(3).specrec,'log')
  
    %hca = irf_panel('fz by largest K');  
    %irf_spectrogram(hca,G1b.tsReduced(3).specrec,'log')
  
    hlinks = linkprop(h([2 4:end]),{'CLim','YLim'});
    
  
    irf_plot_axis_align
    irf_timeaxis(h)
    if doPrint, cn.print(str_print); end
  
    %% Plot, cold components, 1
    figure(202)
    h = irf_plot(6);

    Ks = unique(GM.nComp);
    flog = 'log';
  
    hca = irf_panel('omni fpi');
    irf_spectrogram(hca,GM.userdata.pdist_derived.omni.deflux.specrec)
  
    hca = irf_panel('fz fpi');
    irf_spectrogram(hca,GM.userdata.pdist_derived.fz.specrec,flog)  
  

    hca = irf_panel('fz by max K');  
    irf_spectrogram(hca,GMaxK.tsReduced(3,'vg',GM.userdata.pdist_derived.fz.depend{1}(1,:)).specrec,flog)

    hca = irf_panel('BIC/minBIC');
    irf_spectrogram(hca,GM.tsBIC('normalization','min').specrec( ...
      'p_label','log BIC-min(BIC)','f_label','K'),'log')
    hold(hca,'on')
    irf_plot(hca,GM.tsBICKnee('method','min'),'k')
    irf_plot(hca,GM.tsBICKnee('method','chord'),'r')
    hold(hca,'off')
  
    hca = irf_panel('fz by BIC min');
    specrec = GBICmin.tsReduced(3,'vg',GM.userdata.pdist_derived.fz.depend{1}(1,:)).specrec('p_label',{'PSD m/m^4','K - BIC min'});
    irf_spectrogram(hca,specrec,flog)
  
    hca = irf_panel('fz by BIC chord');
    specrec = GBICchord.tsReduced(3,'vg',GM.userdata.pdist_derived.fz.depend{1}(1,:)).specrec('p_label',{'PSD m/m^4','K - BIC chord'});
    irf_spectrogram(hca,specrec,flog)
  
  
    %hca = irf_panel('fz by largest K');  
    %irf_spectrogram(hca,G1b.tsReduced(3).specrec,'log')
  
    if 0
      hca = irf_panel('nComp cold');
      irf_plot(hca,{GBICmin.tsNComp,GBICchord.tsNComp},'comp')
      hca.YLim = [min(Ks) max(Ks)];    
    end

    hlinks = linkprop(h([2 3 5 6]),{'CLim','YLim'});
    h(5).CLim = prctile(h(5).Children.CData(:),[1 99]);
    
    irf_plot_axis_align
    irf_timeaxis(h)
  
    %% Plot, scatter
    figure(210)
    Ks = 1:4;

    h = setup_subplots(2,numel(Ks),'horizontal');
    isub = 1;
    %Ks = GM.slice(1,':').nComp;
    for iK = 1:numel(Ks)
      K = Ks(iK);
      hca = h(isub); isub = isub + 1;
      G = GM.slice(':',K).sortBy('T');
      tsN = G.tsN('n',GM.userdata.TSeries.n.data);
      tsT = G.tsT('type','iso');
    
      cellN = cellfun(@(x)x.data,tsN,'UniformOutput',false);
      matN = cat(2,cellN{:});
      cellT = cellfun(@(x)x.data,tsT,'UniformOutput',false);
      matT = cat(2,cellT{:});
      
      hca.ColorOrder = pic_colors('matlab');
      scatter(hca,matN,matT,'o')
      
      if 1 % add a fit
        %%
        hca.ColorOrder = pic_colors('matlab');
        hold(hca,'on')
        for ik = 1:K
          x = (matN(:,ik));
          y = (matT(:,ik));
          idx = isfinite(x) & isfinite(y);
          x = x(idx);
          y = y(idx);
          pow = 3;
          [p,S,mu] = polyfit(x,y,pow);
          x_ = linspace(min(x),max(x),100);
          f = polyval(p,x_,S,mu);
          %f = sum((x_').^(0:(pow-1))*(p-1),2);
          %f = p(1)*x_.^(pow-1) + p(2)*x_.^(pow-2) + p(3)*x_.^(pow-3) + p(4)*x_.^(pow-4);
          f(f<min(y)) = NaN;
          f(f>max(y)) = NaN;
          plot(hca,x_,f,'--','linewidth',1)
        end
        hold(hca,'off')
      end

      hca.XLabel.String = 'n (cm^{-3})';
      hca.YLabel.String = 'T (eV)';
    end  

    for iK = 1:numel(Ks)
      K = Ks(iK);
      hca = h(isub); isub = isub + 1;
      G = GM.slice(':',iK).sortBy('T');
      tsV = G.tsV;
  
      cX = cellfun(@(x)x.x.data,tsV,'UniformOutput',false);
      mX = cat(2,cX{:});
      cY = cellfun(@(x)x.z.data,tsV,'UniformOutput',false);
      mY = cat(2,cY{:});
      scatter(hca,mX,mY,'o')
      hca.XLabel.String = 'v_x (km/s)';
      hca.YLabel.String = 'v_z (km/s)';
    end  
    
    c_eval('h(?).Box = ''on'';',1:numel(h))
    c_eval('h(?).FontSize = 12;',1:numel(h))
  
    %% Plot 2, using .byNComp
    h = setup_subplots(2,6,'vertical');
    isub = 1;
  
    for iK = 1:K
      G = GM.slice(':',iK).sortBy('T');
      data = G.byNComp('n',GM.userdata.TSeries.n.data);
  
      hca = h(isub); isub = isub + 1;
      scatter(hca,data.n,data.T,'o')
      %hca.XScale = 'log';
      %hca.YScale = 'log';
  
      hca = h(isub); isub = isub + 1;
      scatter(hca,data.V(:,:,1),data.V(:,:,3),'o')
    end
  
    c_eval('h(?).Box = ''on'';',1:numel(h))
    c_eval('h(?).FontSize = 12;',1:numel(h))
  end
end


