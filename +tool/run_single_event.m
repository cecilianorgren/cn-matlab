% run_single_event
% Specify data folder
%cd /Users/Cecilia/Data/Cluster/20070417/
%savePath = '/Users/Cecilia/Research/LH2/20070417/singles/';
savePath = '/Users/Cecilia/Research/LH2/20080422/singles/ll/';

% Input parameters
sc = 1;
%tint = toepoch([2007 04 17 15 34 22.02; 2007 04 17 15 34 22.10])';
%tint = toepoch([2007 04 17 15 32 50.70; 2007 04 17 15 32 51.00])';
tint_last_cross = toepoch([2008 04 22 17 50 00;2008 04 22 18 15 00])';
tint_ms_msp = toepoch([2008 04 22 17 50 00;2008 04 22 17 57 00])';
tint_msp_ms = toepoch([2008 04 22 17 35 30;2008 04 22 17 37 30])';
tint_zoom = tint_msp_ms;

csys = 'gsm';
flim = 0.5;
v=linspace(20,800,30);

% Obtain correlation parameters
tool.single_event
Vi_loc = irf_resamp(gsmVi1,mean(tint),'nearest'); Vi_loc=Vi_loc(2:4);

doPrint = 0;
if 1
    %%
    
    h=irf_plot(6,'newfigure');
    hca = irf_panel('B'); c_eval('irf_plot(hca,irf_tlim(gsmB?,tint_zoom));',sc); hca.YLabel.String = 'B [nT]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('Vi'); c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom));',sc); hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('E'); c_eval('irf_plot(hca,irf_tlim(gsmE?,tint_zoom));',sc); hca.YLabel.String = 'E [mV/m]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    irf_zoom(h(1:3),'x',tint_zoom)
    irf_pl_mark(h(1:3),tint)
    delete(h(4))
    hca = axes('Position',h(5).Position); delete(h(5)); h(5) = hca;
    irf_match_phibe_vis(hca,'best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc,max(phi_B(:,2))/max(Bz(:,2)));
    
    irf_legend(h(5),['E_n = ' num2str(avEn,'%.f') ' mV/m   E_k = ' num2str(avEk,'%.f') ' mV/m   B_0 = [' num2str(irf_norm(z(1,:)),'%.2f') ']\times' num2str(B0,'%.f') ' nT' ],[0.02 -.5])
    irf_legend(h(5),['v_{ExB} = [' num2str(irf_norm(vExB),'%.2f') ']\times' num2str(irf_abs(vExB,1),'%.f') ' km/s' ],[0.02 -.75])
    irf_legend(h(5),['v-v_{ExB} = [' num2str(irf_norm(direction*velocity-vExB),'%.2f') ']\times' num2str(irf_abs(direction*velocity-vExB,1),'%.f') ' km/s    n = ' num2str(n_loc,'%.1f') ' cc' ],[0.02 -1.0])
    irf_legend(h(5),['\theta_{v>v_{ExB}} = ' num2str(acosd(direction*irf_norm(vExB)'),'%.f') ' ^{\circ}'],[0.02 -1.25])
    %text(h(5).XLim(1),h(5).YLim(2)*1,'a')
    delete(h(6))
    irf_plot_zoomin_lines_between_panels(h(3),h(5))
    title(h(1),['C' num2str(sc) '  ' upper(csys)])
    irf_zoom(h(1:3),'x',tint_zoom)
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_inclov_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
    if 0 % plot small direction plot
        %%
        %figure; 
        qh = subplot(1,1,1);
        qX = [0 0 0]; qU = [x(i_dir,1) y(i_dir,1) z(i_dir,1)];        
        qY = [0 0 0]; qV = [x(i_dir,2) y(i_dir,2) z(i_dir,2)];
        qZ = [0 0 0]; qW = [x(i_dir,3) y(i_dir,3) z(i_dir,3)];
        quiver3(qh,qX,qY,qZ,qU,qV,qW)
        qh.XLabel.String = 'x';
        qh.YLabel.String = 'y';
        qh.ZLabel.String = 'z';
        qh.Title.String = 'GSM';
        
        vratio = avExB/velocity;
        qh.YLim = vratio*[-1 1]; qh.XLim = vratio*[-1 1]; qh.ZLim = vratio*[-1 1];
        axis square
        text(qU(1),qV(1),qW(1),['k, v_{ph}=' num2str(velocity,'%.0f') ' km/s'])
        text(qU(2),qV(2),qW(2),'n')
        text(qU(3),qV(3),qW(3),'B')        
        
        hold(qh,'on')
        quiver3(0,0,0,vratio*avExBdir(1),vratio*avExBdir(2),vratio*avExBdir(3))       
        text(vratio*avExBdir(1),vratio*avExBdir(2),vratio*avExBdir(3),['v_{ExB}=' num2str(avExB,'%.0f') ' km/s'])        
        quiver3(0,0,0,x(i_dir,1)-vratio*avExBdir(1),x(i_dir,2)-vratio*avExBdir(2),x(i_dir,3)-vratio*avExBdir(3))       
        text(x(i_dir,1)-vratio*avExBdir(1),x(i_dir,2)-vratio*avExBdir(2),x(i_dir,3)-vratio*avExBdir(3),['v_{ph}-ExB=' num2str(irf_abs(direction*velocity-vExB,1),'%.0f') ' km/s'])        
        
        hold(qh,'off')
    end
else
    irf_match_phibe_vis('best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc);
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
end
if 0
%% Visualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass);      
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'C' num2str(sc) '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_dir_.gif'],'DelayTime',0.01,'LoopCount',inf);
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath 'C' num2str(sc) '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);


%% Add the event and save data to TimeTable
comment = 'very good!';
toSave = 'single';
tool.add_event
end