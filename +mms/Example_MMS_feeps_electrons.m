
%   Example code to plot FEEPS-electron burst-l2 data
%       - The current burst-l2 data is not helpful to see the pitch angle distributino. A better way is to ask
%         help from FEEPS team, e.g. Drew Turner (drew.lawson.turner@gmail.com).
%       - Refer to FEEPS and EPD documents for better understanding of the instrument.
%       - Wenya, 2016-03-11, irfu.

%%  1. basic
    ic = 1;
    Tint  = tint;%irf.tint('2015-09-08T10:28:34Z/2015-09-08T10:29:50Z');
%%  2. load
%       2.1. load data
    c_eval('feeps_obj = dataobj(''/Volumes/Nexus/data/mms?/feeps/brst/l2/electron/2017/07/11/mms?_feeps_brst_l2_electron_20170711222923_v5.5.1.cdf'');',ic)            % [ion moments]
    eElow = get_variable(feeps_obj, 'electron_energy_lower_bound'); % 1 x 16
    eEupp = get_variable(feeps_obj, 'electron_energy_upper_bound'); % 1 x 16
    eEnergy = get_variable(feeps_obj, 'electron_energy'); eEnergy = eEnergy.data; % 1 x 16
    
    c_eval('eSpin = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_spin'');', ic); % nt x 1, 'The number of spacecraft rotations'
    c_eval('eSpinIntegrationSectors = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_integration_sectors'');', ic); % nt x 1, all is 1
    c_eval('eSpinDuration = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_spin_duration'');', ic); % nt x 1, close to 19.6
    c_eval('ePitchAngle = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_pitch_angle'');', ic); % nt x 18 , range [0 180], pitch_angle of different sectors
    c_eval('eSpinSectNum = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_spinsectnum'');', ic); % nt x 1, range [0 63]
   
    
    % Total number of sectors is 18
    c_eval('eTit1 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_1'');', ic); % nt x 16 (nE)
    c_eval('eTit2 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_2'');', ic);
    c_eval('eTit3 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_3'');', ic);
    c_eval('eTit4 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_4'');', ic);
    c_eval('eTit5 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_5'');', ic);
    c_eval('eTit9 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_9'');', ic);
    c_eval('eTit10 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_10'');', ic);
    c_eval('eTit11 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_11'');', ic);
    c_eval('eTit12 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_12'');', ic);
    c_eval('eBit1 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_1'');', ic);
    c_eval('eBit2 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_2'');', ic);
    c_eval('eBit3 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_3'');', ic);
    c_eval('eBit4 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_4'');', ic);
    c_eval('eBit5 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_5'');', ic);
    c_eval('eBit9 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_9'');', ic);
    c_eval('eBit10 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_10'');', ic);
    c_eval('eBit11 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_11'');', ic);
    c_eval('eBit12 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_12'');', ic);
    c_eval('ePA = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_pitch_angle'');', ic);
%       2.2. tlim data
    c_eval('eTit? = eTit?.tlim(Tint);', [1:5, 9:12]);   
    c_eval('eBit? = eBit?.tlim(Tint);', [1:5, 9:12]);
    ePitchAngle = ePitchAngle.tlim(Tint);
    ePA = ePA.tlim(Tint);
    %%
    nT = eTit1.length;
    nEnergy = numel(eEnergy); % 16
    nPitchangles = size(ePitchAngle.data,2); % 18
    eFeepsData = nan(nT,nEnergy,nPitchangles);
    eFeepsData(:,:,1) = eTit1.data;
    eFeepsData(:,:,2) = eTit2.data;
    eFeepsData(:,:,3) = eTit3.data;
    eFeepsData(:,:,4) = eTit4.data;
    eFeepsData(:,:,5) = eTit5.data;
    eFeepsData(:,:,6) = eTit9.data;
    eFeepsData(:,:,7) = eTit10.data;
    eFeepsData(:,:,8) = eTit11.data;
    eFeepsData(:,:,9) = eTit12.data;
    eFeepsData(:,:,10) = eBit1.data;
    eFeepsData(:,:,11) = eBit2.data;
    eFeepsData(:,:,12) = eBit3.data;
    eFeepsData(:,:,13) = eBit4.data;
    eFeepsData(:,:,14) = eBit5.data;
    eFeepsData(:,:,15) = eBit9.data;
    eFeepsData(:,:,16) = eBit10.data;
    eFeepsData(:,:,17) = eBit11.data;
    eFeepsData(:,:,18) = eBit12.data;
    eFeepsData(eFeepsData==0) = NaN;
    eEnergyTable = repmat(torow(eEnergy),nT,1);
    [sortedPitchangles,iSortedPitchangles] = sort(ePitchAngle.data,2);
    eFeepsDataSorted = eFeepsData(iSortedPitchangles);
    c_eval('eFeepsPitch? = irf.ts_pitchangles(eTit?.time,eFeepsData,eEnergyTable,ePitchAngle.data)',ic)
    c_eval('eFeepsPitch? = irf.ts_pitchangles(eTit?.time,eFeepsData,eEnergyTable,ePitchAngle.data)',ic)
    
    
%%  3. make
%       3.1. energy bins
    eenergy = (eElow.data + eEupp.data)/2.;
    
%       3.2. structure for Top
    c_eval('speTit? = struct(''t'', eTit?.time.epochUnix);', [1:5, 9:12]);
    c_eval('speTit?.p = double(eTit?.data);', [1:5, 9:12]);
    c_eval('speTit?.p_label = {''intensity''};', [1:5, 9:12]);
    c_eval('speTit?.f_label = {''Energy''};', [1:5, 9:12]);
    c_eval('speTit?.f = double(eenergy);', [1:5, 9:12]);
    
%       3.3. structure for Bottom
    c_eval('speBit? = struct(''t'', eBit?.time.epochUnix);', [1:5, 9:12]);
    c_eval('speBit?.p = double(eBit?.data);', [1:5, 9:12]);
    c_eval('speBit?.p_label = {''intensity''};', [1:5, 9:12]);
    c_eval('speBit?.f_label = {''Energy''};', [1:5, 9:12]);
    c_eval('speBit?.f = double(eenergy);', [1:5, 9:12]);

%%  4. plot
%       4.0.
    h = irf_plot(10,'newfigure');
%       4.1. electron Top 5
    hca=irf_panel('Top1');
    irf_spectrogram(hca, speTit1, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H1U [keV]','Interpreter','tex');

%       4.2. electron Top 3
    hca=irf_panel('Top3');
    irf_spectrogram(hca, speTit3, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H2U','Interpreter','tex');
    
%       4.3. electron Top 9
    hca=irf_panel('Top9');
    irf_spectrogram(hca, speTit9, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H5U','Interpreter','tex');

%       4.4. electron Top 11
    hca=irf_panel('Top11');
    irf_spectrogram(hca, speTit11, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H6U','Interpreter','tex'); 
    
%       4.5. electron Top 5
    hca=irf_panel('Top5');
    irf_spectrogram(hca, speTit5, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3T','Interpreter','tex');
    
%       4.6. electron Bottom 5
    hca=irf_panel('Bottom5');
    irf_spectrogram(hca, speBit5, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3B','Interpreter','tex');      
    
%       4.7. electron Bottom 2
    hca=irf_panel('Top2');
    irf_spectrogram(hca, speTit2, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H1D','Interpreter','tex');

%       4.8. electron Top 4
    hca=irf_panel('Top4');
    irf_spectrogram(hca, speTit4, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H2D','Interpreter','tex');     

%       4.9. electron Top 10
    hca=irf_panel('Top10');
    irf_spectrogram(hca, speTit10, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3D','Interpreter','tex');

%       4.10. electron Top 12
    hca=irf_panel('Top12');
    irf_spectrogram(hca, speTit12, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H4D','Interpreter','tex');     
    
%       4.X.
    irf_zoom(h,'y', [20 200])
    irf_zoom(h,'x',Tint)
    irf_plot_axis_align(h)
    %set(h ,'position',[0.15 0.2 0.66 0.71]);    
    %irf_colormap('poynting')    
    %add_position(h(end),Rgse)
    %xlabel(h(end),'');
    tmpts0 = Tint.start.utc;
    tmpts1 = Tint.stop.utc;
    title(h(1),['MMS' int2str(ic) ' ' tmpts0(12:22) '~' tmpts1(12:22) ]);

%%