%% Example code to perform Walen test; only for burst mode MMS data.
%   Reference: Phan et al. [AG, 2004]

%%  1. basic
    ic = 1;
    Jsign = 1;         % 1/-1 for jet direction;    
    time = irf_time('2015-11-12T06:43:00.000Z','utc>epochtt');
    vec1 = [1, 0, 0];       vec2 = [0, 1, 0];       vec3 = [0, 0, 1];       % in GSE
    tint = irf.tint('2015-11-12T07:19:00.000000Z/2015-11-12T07:22:00.000000Z');             % plot
    tint_ref = irf.tint('2015-11-12T07:21:20.000000Z/2015-11-12T07:21:30.000000Z');         % reference region     
    tint_walen = irf.tint('2015-11-12T07:19:00.000000Z/2015-11-12T07:21:10.000000Z');       % test region
    unit_v2w = 1e6 * 0.5 * 1.673e-27 / 1.6e-19;
    
%%  1. basic, different time 
    ic = 1;
    Jsign = -1;         % 1/-1 for jet direction in time!;    
    time = irf_time('2015-11-12T06:40:00.000Z','utc>epochtt');
    vec1 = [0.3229   -0.0940    0.9418];       vec2 = [-0.1780    0.9713    0.1580];       vec3 = [-0.9296   -0.2187    0.2969];       % in GSE         
    tint = irf.tint('2015-11-12T06:40:00.000000Z/2015-11-12T06:45:00.000000Z');             % plot
    tint_ref = irf.tint('2015-11-12T06:44:00.000000Z/2015-11-12T06:44:30.000000Z');         % reference region     
    tint_walen = irf.tint('2015-11-12T06:40:55.000000Z/2015-11-12T06:43:00.000000Z');       % test region
    unit_v2w = 1e6 * 0.5 * 1.673e-27 / 1.6e-19;
    
%%  2. load data
    % 2.1. PSD
    c_eval('ifile? = mms.get_filepath([''mms?_fpi_fast_l2_dis-dist''], time);', ic);
    c_eval('[iPDist?, iPDistError?] = mms.make_pdist(ifile?);', ic); 
    % 2.2. moment
    c_eval('Tint = irf.tint(iPDist?.time.start, iPDist?.time.stop);', ic);
    c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',Tint, ?);',ic);
    c_eval('Vi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint,?);',ic);
    c_eval('Pi? = mms.get_data(''Pi_dbcs_fpi_brst_l2'',Tint,?);',ic);   c_eval('Pi?.units = ''nPa'';', ic);
    c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',Tint, ?);',ic);
    % 2.3. fields
    c_eval('gseB?=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);
    % 2.4. Load defatt files
    c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',ic);
    c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',ic);

%%  3. compute         
    % 3.2. alpha: pressure anisotropy factor
    c_eval('alpha = irf_pres_anis(facPi?fast, gseB?srvy.resample(facPi?fast));', ic);
    
    % 3.3. gse to new123
    c_eval('B123 = irf_newxyz(gseB?srvy, vec1, vec2, vec3);', ic);
    c_eval('Vi123 = irf_newxyz(gseVi?fast, vec1, vec2, vec3);', ic);
    
    % 3.4. reference (MSH) region; in New frame (123);
    Bref = B123.tlim(tint_ref);                     Bref = irf.nanmean(Bref.data, 1);
    Viref = Vi123.tlim(tint_ref);                   Viref = irf.nanmean(Viref.data, 1);
    c_eval('Niref = ni?fast.tlim(tint_ref);', ic);      Niref = irf.nanmean(Niref.data, 1);
    alpharef = alpha.tlim(tint_ref);                alpharef = irf.nanmean(alpharef.data, 1);
    
    % 3.5. Vipred1: delta_B / sqrt(rho1)
    c_eval('B123 = B123.resample(ni?fast);', ic);
    c_eval('Vi123 = Vi123.resample(ni?fast);', ic);
    tmp1 = (B123 - Bref) * 21.8 / sqrt(Niref);
    Vipred1 = tmp1.resample(Vi123) * Jsign + Viref;
    
    % 3.6. Vipred2: B_2 / sqrt(rho2) - B_1 / sqrt(rho1)
    %c_eval('tmp2 = 21.8 * B123.data ./ sqrt([ni?.data ni?.data ni?.data]);', ic);
    %tmp2 = irf.ts_vec_xyz(B123.time, tmp2);
    %Vipred2 = (tmp2 - 21.8 * Bref / sqrt(Niref)) * Jsign + Viref;
    c_eval('tmp2 = 21.8 * (1-alpha) * B123.resample(alpha) / sqrt(Niref) / sqrt(1-alpharef);', ic);       % Phan et al. [AG, 2004]
    Vipred2 = (tmp2 - 21.8 * sqrt(1-alpharef) * Bref / sqrt(Niref)) * Jsign + Viref;
    
    % 3.7. Vipred2: sqrt(1-alpha_2) * B_2 / sqrt(rho2) - sqrt(1-alpha_2) * B_1 / sqrt(rho1)  
    c_eval('tmp3 = 21.8 * sqrt(1-alpha) * B123.resample(alpha) / sqrt(ni?fast);', ic);
    Vipred3 = (tmp3 - sqrt(1 - alpharef) * 21.8 * Bref / sqrt(Niref)) * Jsign + Viref;

    % 3.8. slope & CC
    Vi123w = Vi123.tlim(tint_walen);             Vipredw1 = Vipred1.tlim(tint_walen);    
    Vipredw2 = Vipred2.tlim(tint_walen);    Vipredw3 = Vipred3.tlim(tint_walen);    
    c_eval('p? = polyfit(Vipredw2.?.data, Vi123w.?.data, 1);     slope2? = p?(1);', ['x', 'y', 'z']);
    c_eval('corr? = corrcoef(Vipredw2.?.data, Vi123w.?.data);    cc2? = corr?(1, 2);', ['x', 'y', 'z']);    
   
%%  4. plot
    % 4.0. figure setting
    h = irf_plot(7,'newfigure');
    xSize=800; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.129;
    c_eval('set(h(?), ''position'', [0.12 0.960-? * ywidth xwidth ywidth]);', 1: 7);
    
    % 4.1. Bxyz
    h(1) = irf_panel('Bgse');
    c_eval('irf_plot(h(1), gseB?srvy, ''LineWidth'', 1.5);', ic);
    ylabel(h(1),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{x} ',' B_{y} ',' B_{z}'},[0.88 0.15])
    irf_legend(h(1),'(a)',[0.99 0.98],'color','k')
    
    % 4.2. Ni & Ne
    h(2) = irf_panel('ni');
    c_eval('irf_plot(h(2), ni?fast, ''LineWidth'', 1.5);', ic);
    hold(h(2), 'on');
    c_eval('irf_plot(h(2), ne?fast, ''LineWidth'', 1.5);', ic);
    hold(h(2), 'off');
    ylabel(h(2),{'N','[cm^{-3}]'},'Interpreter','tex');
    irf_legend(h(2),{'N_{i} ',' N_{e} '},[0.88 0.15])
    irf_legend(h(2),'(b)',[0.99 0.98],'color','k')    
        
    % 4.3. ion energy flux
    h(3) = irf_panel('iEnflux');
    c_eval('irf_spectrogram(h(3), iPDist?fast.deflux.omni.specrec, ''log'');', ic);
    h(3).YScale = 'log';
    h(3).YTick = 10.^[1 2 3 4];
    irf_legend(h(3),'(c)',[0.99 0.98],'color','k')
    ylabel(h(3),{'W_i','[eV]'},'Interpreter','tex');

    % 4.4. B123
    h(4) = irf_panel('B123');
    irf_plot(h(4), B123, 'LineWidth', 1.5);
    ylabel(h(4),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(4),{'B_{1} ',' B_{2} ',' B_{3}'},[0.88 0.15])
    irf_legend(h(4),'(d)',[0.99 0.98],'color','k')
    % vectors legend
    vec1_str = ['[', num2str(vec1(1),'% .2f'), ',', num2str(vec1(2),'% .2f'), ',', num2str(vec1(3),'% .2f'), ']'];
    vec2_str = ['[', num2str(vec2(1),'% .2f'), ',', num2str(vec2(2),'% .2f'), ',', num2str(vec2(3),'% .2f'), ']'];
    vec3_str = ['[', num2str(vec3(1),'% .2f'), ',', num2str(vec3(2),'% .2f'), ',', num2str(vec3(3),'% .2f'), ']'];
    irf_legend(h(4), vec1_str, [1.003 0.8],'color','k', 'FontSize', 11)
    irf_legend(h(4), vec2_str, [1.003 0.5],'color','b', 'FontSize', 11)
    irf_legend(h(4), vec3_str, [1.003 0.1],'color','r', 'FontSize', 11)
       
    % 4.5. Vi123.x
    h(5) = irf_panel('Vi1');
    irf_plot(h(5), Vipredw2.x, 'LineWidth', 1.5, 'color', 'r');
    hold(h(5), 'on');
    irf_plot(h(5), Vi123.x, 'LineWidth', 1.5, 'color', 'k');
    irf_plot(h(5), Vipred2.x, '--', 'LineWidth', 1.5, 'color', 'r');    
    hold(h(5), 'off');
    ylabel(h(5),{'V_1','[km/s]'},'Interpreter','tex');
    irf_legend(h(5),{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster');
    irf_legend(h(5),'(e)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2x, '% .2f')];
    ccx_str = ['cc=', num2str(cc2x, '% .2f')];
    irf_legend(h(5), slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(h(5), ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.6. Vi123.y
    h(6) = irf_panel('Vi2');
    irf_plot(h(6), Vipredw2.y, 'LineWidth', 1.5, 'color', 'r');
    hold(h(6), 'on');
    irf_plot(h(6), Vi123.y, 'LineWidth', 1.5, 'color', 'k');
    irf_plot(h(6), Vipred2.y, '--', 'LineWidth', 1.5, 'color', 'r');    
    hold(h(6), 'off');
    ylabel(h(6),{'V_2','[km/s]'},'Interpreter','tex');
    irf_legend(h(6),{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster')
    irf_legend(h(6),'(f)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2y, '% .2f')];
    ccx_str = ['cc=', num2str(cc2y, '% .2f')];
    irf_legend(h(6), slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(h(6), ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.5. Vi123.z
    h(7) = irf_panel('Vi3');
    irf_plot(h(7), Vipredw2.z, 'LineWidth', 1.5, 'color', 'r');
    hold(h(7), 'on');
    irf_plot(h(7), Vi123.z, 'LineWidth', 1.5, 'color', 'k');
    irf_plot(h(7), Vipred2.z, '--', 'LineWidth', 1.5, 'color', 'r');    
    hold(h(7), 'off');
    ylabel(h(7),{'V_3','[km/s]'},'Interpreter','tex');
    irf_legend(h(7),{'fpi ',' pred '},[0.88 0.15], 'color', 'cluster')
    irf_legend(h(7),'(g)',[0.99 0.98],'color','k')     
    % slope & cc
    slopex_str = ['slope=', num2str(slope2z, '% .2f')];
    ccx_str = ['cc=', num2str(cc2z, '% .2f')];
    irf_legend(h(7), slopex_str, [1.015 0.65],'color','k', 'FontSize', 12);
    irf_legend(h(7), ccx_str, [1.015 0.25],'color','k', 'FontSize', 12);
    
    % 4.X. global
    load('caa/cmap.mat');
    colormap(h(3),cmap); 
    title(h(1), strcat('MMS-', num2str(ic)));
    irf_plot_axis_align(h(1 : 7))
    irf_zoom(h(1 : 7),'x', tint);
    irf_pl_mark(h(1: 7), tint_ref, 'red')
    irf_pl_mark(h(5: 7), tint_walen, 'yellow'); %, 'LineStyle', '--', 'LineWidth', 1.5)    

%%

