% Figure for for Maria Elenas paper
cd /Users/Cecilia/Data/Cluster/20070902/
savePath = '/Users/Cecilia/Research/LH2/20070902/';

% Time intervals
tint = toepoch([2007 09 02 15 47 43.2; 2007 09 02 15 47 49.9])';
tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';
tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 36.9])';
em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])';
es_tint = toepoch([2008 04 22 17 36 45.2;2008 04 22 17 36 46.1])';
tint_out = toepoch([2008 04 22 17 35 30;2008 04 22 17 37 30])';
tint = es_tint;% Correlation
sc = 1;
csys = 'gsm';

c_eval('R_loc = irf_resamp(gsmR?,tint(1),''nearest'');',sc); R_loc = R_loc(2:4);
c_eval('Vi_loc = irf_resamp(gsmVi?,tint(1),''nearest'');',sc); Vi_loc = Vi_loc(2:4);
c_eval('Ti_loc = irf_resamp(Ti?,tint(1),''nearest'');',sc); Ti_loc = Ti_loc(2);
c_eval('Te_loc = irf_resamp(parTe?,tint(1),''nearest'');',sc); Te_loc = Te_loc(2);


% phi colors
colors_phi = [0.7500    0.1250    0.0980; 
              0.9290    0.6940    0.1250
              0.8500    0.3250    0.0980; 
              0 0 1; 
              1    0.3250    0.0980;   
              0.9290    0.6940    0.1250
              0.8500    0.3250    0.0980; 
              0.9290    0.6940    0.1250];
% length colors              
lengthcolors = [0    0.4470    0.7410; 
                0.4660    0.6740    0.1880; 
                0.4940    0.1840    0.5560;
                0.3010    0.7450    0.9330;
                0.8500    0.3250    0.0980; 
                0.9290    0.6940    0.1250; 
                1 0 0];
                
doPrint = 0;
%% Waveform figure
if 1 % Make main figure
    %tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 34])';
    %tint_out = tint + 1*[-60 60];
    tint_zoom = tint_out;
    h = irf_plot(6,'newfigure');
    % B
    hca = irf_panel('B'); c_eval('irf_plot(hca,irf_tlim(irf_abs(gsmB?),tint_zoom+[-4 4]));',sc); 
    hca.YLabel.String = 'B [nT]'; irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.99 0.95])
    grid(hca,'off'); 
    %hca.YLim = [0 50];
    hca.Title.String = 'Cluster 3';
    %if hca.YLim(1)>0; hca.YLim = hca.YLim.*[0 1]; end; hca.YTick = [10:10:70];
    irf_legend(hca,'a)',[0.02 0.95]);
    
    % n
    %hca = irf_panel('n'); c_eval('irf_plot(hca,{irf_tlim(peaNe?,tint_zoom+[-8 8]),irf_tlim([scpNe?(:,1) scpNe?(:,2)*0.1],tint_zoom+[-4 4])},''comp'');',sc); 
    %hca.YLabel.String = 'n_e [cm^{-3}]'; irf_legend(hca,{'PEACE','EFW'},[0.99 0.95])
    %grid(hca,'off'); 
    %hca.YTick = [0:0.2:10];
    %hca.YLim = [0 50];
    
    % Vi
    hca = irf_panel('Vi'); c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom+[-4 4]));',sc); 
    %hca = irf_panel('Vi'); c_eval('irf_plot(hca,{irf_tlim(gsmExB?,tint_zoom+[-4 4]),irf_tlim(gsmVi?,tint_zoom+[-4 4])},''comp'');',sc); 
    hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'V_x','V_y','V_z'},[0.99 0.95])
    grid(hca,'off')
    irf_legend(hca,'b)',[0.02 0.95]);
    % E
    hca = irf_panel('E'); c_eval('irf_plot(hca,irf_tlim(gsmE?,tint_zoom+[-4 4]));',sc); 
    hca.YLabel.String = 'E [mv/m]'; irf_legend(hca,{'E_x','E_y','E_z'},[0.99 0.95])
    grid(hca,'off')
    hca.YLim = [-50 50];
    hca.YTick = [-40:20:40]; 
    irf_legend(hca,'c)',[0.02 0.95]);
    
    irf_zoom(h(1:4),'x',tint_zoom)
    irf_pl_mark(h(1:4),tint)
    hca = irf_panel('delete for space'); %delete(hca)    
    
    grid(hca,'off')    
    irf_timeaxis(h(3))
    %colors_phi = [0 0 0; 1 0 0];
    
    phi_labels_all = {{'\Phi_{\delta E_{k}}','\Phi_{\delta B_{||}}'},...
                      {'\Phi_{\delta E{k}}','\Phi_{\delta B{||}}'},...
                      {'\Phi_{\delta E_{\perp}}','\Phi_{\delta B_{||}}'},...
                      {'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},...
                      {'\Phi_{\delta E}','\Phi_{\delta B}'}};
    phi_labels = phi_labels_all{5};
                    
    
    % low f phi
    flim = 0.1;
    tool.single_event       
    hca = irf_panel('Phi: low f');
    hca.ColorOrder = colors_phi;
    % move panel slightly upwards to make room for km labels on panel below
    hca.Position = hca.Position + [0 0.02 0 0];
    irf_plot(hca,{phi_E(:,[1 1+i_v]),phi_B},'comp')    
    hca.YLabel.String = '\Phi [V]'; 
    %hleg = irf_legend(hca,{'\Phi_E','\Phi_B'},[0.99 0.95]);
    hleg = irf_legend(hca,phi_labels,[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E}','\Phi_{\delta B}'},[0.99 0.95]);
    hleg(1).Color = colors_phi(1,:);
    hleg(2).Color = colors_phi(2,:);
    irf_legend(hca,{['     v_{ph}=' num2str(velocity,'%.0f') ' km/s']},[0.02 0.95])
    irf_legend(hca,{['\omega > ' num2str(flim,'%.2f') ' \Omega_{LH}']},[0.02 0.1])
    irf_zoom(hca,'x',tint)
    %irf_zoom(hca,'y',[-25 25])    
    axtop = add_length_on_top(hca,velocity,1);
    axtop.XLabel.String = 'Length [km]';    
    % add B labels on right
    yratio = phi_B(100,2)/Bz(100,2);            
    axtop.YLim = hca.YLim/yratio;
    axtop.YLabel.String = '\delta B_{||} [nT]'; 
    axtop.YTickLabelMode = 'auto';
    axtop.YTickMode = 'auto';
    % mark different physical lengths
    if 0
    T1 = axtop.XTick(3)+irf_plot_start_epoch; mean(tint)-diff(tint)*0.2;
    T0 = axtop.XTick(3)+irf_plot_start_epoch; mean(tint);
    T1 = axtop.XTick(1)+irf_plot_start_epoch;
    T_di = di_loc/velocity; t_di = [T1+[0 T_di]' hca.YLim(2)*0.8*[1 1]']; th = text(t_di(2,1)-irf_plot_start_epoch,t_di(2,2)/yratio*0.85,' d_i','color',lengthcolors(1,:),'horizontalalignment','left');
    T_ri = ri_loc/velocity; t_ri = [T1+[0 T_ri]',hca.YLim(2)*0.9*[1 1]']; th = text(t_ri(2,1)-irf_plot_start_epoch,t_ri(2,2)/yratio*0.9,' r_i','color',lengthcolors(2,:),'horizontalalignment','left');
    T_reri = sqrt(ri_loc*re_loc)/velocity; t_reri = [T1+[0 T_reri]',hca.YLim(2)*0.7*[1 1]']; th = text(t_reri(2,1)-irf_plot_start_epoch+0,t_reri(2,2)/yratio*0.75,' (r_er_i)^{1/2}','color',lengthcolors(3,:),'horizontalalignment','left' );
    %T_re = re_loc/velocity; t_re = [T1+[0 T_re]',hca.YLim(2)*0.6*[1 1]']; th = text(t_re(2,1)-irf_plot_start_epoch+0,t_re(2,2)/yratio,'\rho_e','color',lengthcolors(4,:));
    
    hold(hca,'on')
    
    irf_plot(hca,t_di,'linewidth',2,'color',lengthcolors(1,:))
    irf_plot(hca,t_ri,'linewidth',2,'color',lengthcolors(2,:))
    irf_plot(hca,t_reri,'linewidth',2,'color',lengthcolors(3,:))
    %irf_plot(hca,t_re,'linewidth',2)
    hold(hca,'off')  
    end
    grid(hca,'off')    
    hca.XTickLabel = [];
    hca.XLabel.String = '';  
    irf_legend(hca,'d)',[0.02 0.95]);
    velocity;
    
    % Compare Ek and Bz amplitudes
    maxEk = max(dEk(:,2));
    maxEn = max(dEn(:,2));
    maxBz = max(Bz(:,2));
    disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
    
    % high f phi
    flim = 0.5;
    tool.single_event   
    hca = irf_panel('Phi: high f');
    hca.ColorOrder = colors_phi;
    % move panel slightly downwards to make room for km labels on top
    hca.Position = hca.Position + [0 -0.02 0 0];
    irf_plot(hca,{phi_E(:,[1 1+i_v]),phi_B},'comp')
    hca.YLabel.String = '\Phi [V]'; 
    hleg = irf_legend(hca,phi_labels,[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E}','\Phi_{\delta B}'},[0.99 0.95]);
    hleg(1).Color = colors_phi(1,:);
    hleg(2).Color = colors_phi(2,:);
    irf_legend(hca,{['      v_{ph}=' num2str(velocity,'%.0f') ' km/s']},[0.02 0.95])
    irf_legend(hca,{['\omega > ' num2str(flim,'%.2f') ' \Omega_{LH}']},[0.02 0.1])
    irf_zoom(hca,'x',tint)
    %irf_zoom(hca,'y',[-15 15]) 
    axtop = add_length_on_top(hca,velocity,0.25);    
    % add B labels on right
    yratio = phi_B(100,2)/Bz(100,2);            
    axtop.YLim = hca.YLim/yratio;
    axtop.YLabel.String = '\delta B_{||} [nT]'; 
    axtop.YTickLabelMode = 'auto';
    axtop.YTickMode = 'auto';
    % mark different physical lengths
    hca = irf_panel('Phi: high f');
    %T1 = axtop.XTick(4)+irf_plot_start_epoch; mean(tint)-diff(tint)*0.2;
    %T0 = axtop.XTick(4)+irf_plot_start_epoch; mean(tint);
    %2
    %T_di = di_loc/velocity; t_di = [T1+[0 T_di]' hca.YLim(2)*0.8*[1 1]']; %th = text(t_di(2,1)-irf_plot_start_epoch,t_di(2,2)/yratio,'d_i')
    %T_re = re_loc/velocity; t_re = [T1+[0 T_re]',hca.YLim(2)*0.9*[1 1]'];
    %T_reri = sqrt(ri_loc*re_loc)/velocity; t_reri = [T1+[0 T_reri]',hca.YLim(2)*0.7*[1 1]'];
    %T_ri = ri_loc/velocity; t_ri = [T1+[0 T_ri]',hca.YLim(2)*0.9*[1 1]'];
    %hold(hca,'on')
    %irf_plot(hca,t_ri,'linewidth',2)
    %irf_plot(hca,t_di,'linewidth',2)
    %irf_plot(hca,t_reri,'linewidth',2)
    %irf_plot(hca,t_re,'linewidth',2,'color',lengthcolors(4,:)); t_re = [T1+[0 T_re]',hca.YLim(2)*0.9*[1 1]']; th = text(t_re(2,1)-irf_plot_start_epoch+0,t_re(2,2)/yratio*0.92,' r_e','color',lengthcolors(4,:));
    %hold(hca,'off')    
    grid(hca,'off')
    irf_legend(hca,'e)',[0.02 0.95]);
    velocity;
    
    irf_plot_axis_align
    irf_plot_zoomin_lines_between_panels(h(3),h(5))
    hca = irf_panel('delete for space'); delete(hca) 
    
    % Compare Ek and Bz amplitudes
     maxEk = max(dEk(:,2));
    maxEn = max(dEn(:,2));
    maxBz = max(Bz(:,2));
    disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
end