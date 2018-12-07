tint = irf.tint('2015-11-12T07:19:05.00Z/2015-11-12T07:19:44.00Z'); %20151112071854
ic = 1;
localuser = datastore('local','user');

% Load datastore
mms.db_init('local_file_db','/Volumes/Nexus/data/');
db_info = datastore('mms_db');   

% Make event directory
eventPath = ['/Users/' localuser '/GoogleDrive/Research/Outreach/'];

%% Load data
units = irf_units;
% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);

% Particle moments
% Density
disp('Loading density...');
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);

% Velocity
disp('Loading bulk velocities...');
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)

% Derived quantities
c_eval('gseJ? = units.e*ne?*1e6*(gseVi?.resample(ne?)-gseVe?)*1e3*1e9;',ic)
c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,12);',ic)

%% Make plot
doPrint = 0;
doGif = 1;

t_cluster = 4;
t_mms = 0.03;
t_ratio = t_cluster/t_mms;

npanels = 4;
%h = irf_plot(npanels);
i_frame = 0;
t_step = 3;
%%
for it = 1;%t_ratio:-t_step:1%1:t_step:t_ratio
  i_frame = i_frame + 1;
  if it == 1
    timeline = ePDist1.time;
  else
    timeline = ePDist1.time(1):(t_mms*it):ePDist1.time(end);
  end
  %timeline = tint(1):(t_mms*it):tint(2);
  if ~isequal(ePitch1.time,timeline)
    ePitch = ePitch1.resample(timeline);
    depend1 = resample(irf.ts_scalar(ePitch1.time,ePitch.depend{1}),timeline);
    ePitch.depend{1} = depend1.data;
  else
    ePitch = ePitch1;
  end
  if ~isequal(ePDist1.time,timeline)
    ePDist = ePDist1.resample(timeline);
    depend1 = resample(irf.ts_scalar(ePDist1.time,ePDist.depend{1}),timeline);
    ePDist.depend{1} = depend1.data;
  else
    ePDist = ePDist1;
  end

  if 1
    hca = irf_panel('B');
    irf_plot(hca,gseB1);
    hca.YLabel.String = {'Magnetic','field', '(nT)'};
  end
  if 1
    hca = irf_panel('J');
    irf_plot(hca,gseJ1.resample(timeline));
    hca.YLabel.String = {'Current','density', '(nA/m^2)'};
    hca.YLabel.Interpreter = 'tex';
    hca.YLim = [-600 900];
  end
  if 1
    hca = irf_panel('electron omni def');  
    [~,hcb1] = irf_spectrogram(hca,ePDist.deflux.omni.specrec);
    hca.YLabel.String = {'Electron','energy', '(eV)'};
    hca.YLabel.Interpreter = 'tex';
    hca.YScale = 'log';
    hca.CLim = [3 8.3];
    %hca.CLim = [-35 -27];
    hca.YTick = 10.^[2:10];
    hca.YLim = [20 30000];
  end
  if 1
    hca = irf_panel('electron pitch def');  
    [~,hcb2] = irf_spectrogram(hca,ePitch.deflux.specrec('pa'));
    hca.YLabel.String = {'Electron','pitch angle','(deg.)'};
    hca.YLabel.Interpreter = 'tex';    
    hca.CLim = [6.8 7.9];
  end

      
  irf_zoom(h,'x',tint+[1 -1])
  irf_plot_axis_align
  
  for ip = 1:npanels
    h(ip).Position(3) = 0.7;
  end
  hcb1.Position(1) = 0.88;
  hcb2.Position(1) = 0.88;
  
  if doPrint
    cn.print(sprintf('downsampling_mp_dt=%g',it*t_mms),'path',eventPath)
  end
  if doGif
    drawnow;    
    f = getframe(gcf);
    A(i_frame) = f;
    if it == 1 % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,20) = 0;
    else
        im(:,:,1,i_frame) = rgb2ind(f.cdata,map,'nodither');
    end
  end  
end

if doGif
  %movie(fg,A,20);
  imwrite(im,map,sprintf('%sgifs/downsampling_t_step%g_reverse_1.gif',eventPath,t_step),'DelayTime',0.1,'LoopCount',0)
end
