% small fluxrope?
savePath = '/Users/Cecilia/Research/LH2/20070417/small_fluxrope/';
%% Overview plot
sc=1:4;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
tintlim = toepoch([2007 04 17 15 32 00;2007 04 17 15 33 00])';
h = irf_plot(5,'newfigure');
isub = 1;

if 1 % B3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB3);
    hca.YLabel.String = 'B3 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB4);
    hca.YLabel.String = 'B4 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE3);
    hca.YLabel.String = 'E3 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE4);
    hca.YLabel.String = 'E4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % Vi3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3);
    hca.YLabel.String = 'Vi4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end

title(h(1),'')
irf_zoom(h,'x',tintlim)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];

%% Moving average
ffiltlim = 1;
nTlh = 3; 
vrange = [50 800];
sc = 3;
tintlim = toepoch([2007 04 17 15 32 18;2007 04 17 15 32 24])'; % C3#1
%tintlim = toepoch([2007 04 17 15 32 38;2007 04 17 15 32 44])'; % C3#2
%tintlim = toepoch([2007 04 17 15 32 26;2007 04 17 15 32 38])'; % C4
%tintlim = toepoch([2007 04 17 15 32 16;2007 04 17 15 32 44])'; % C3, with a wavegap inbetween
c_eval('[mB,cmaxs,mE,mintEdt,ks,ns,Bs,B0s,cs,mva_v1,mva_v2,mva_v3,mva_l123,EB,tints,ffilts,vs] = tool.moving_average(irf_tlim(gsmE?,tintlim),irf_tlim(gsmB?,tintlim),irf_tlim(flh?,tintlim),ffiltlim,nTlh,peaNe?hf,vrange);',sc);

%% Plot moving average
h = irf_plot(9,'newfigure');
isub = 1;


if 1 % E
    hca = h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_tlim(gsmE?,tintlim));',sc)
    hca.YLabel.String = irf_ssub('E? [mV/m]',sc);
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E maximum variance direction
    hca = h(isub); isub=isub+1;
    irf_plot(hca,mva_v1)    
    hca.YLabel.String = 'Maximum \newline variance \newline direction';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % k from matching
    hca = h(isub); isub=isub+1;
    irf_plot(hca,ks)  
    hca.YLabel.String = 'k';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;
    irf_plot(hca,cmaxs)  
    hca.YLabel.String = 'C_{max}';      
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;   
    specrec.t = tints(:,1);
    specrec.p = cs;
    specrec.f = 1:3:360';
    irf_spectrogram(hca,specrec) 
    hca.YLabel.String = '\theta';
    hca.CLim = [-1 1];
    hc = colorbar('peer',hca);
    hc.YLim = [-1 1];
    hc.YLabel.String = 'C';
    cmap = irf_colormap('poynting');
    colormap(hca,cmap)
end
if 1 % max correlation 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,mE)    
    hca.YLabel.String = '<|\delta E|>';    
end

if 1 % frequencies 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,ffilts); hold(hca,'on');
    c_eval('flhplot = irf_tlim(flh?,tintlim);',sc)
    irf_plot(hca,flhplot); hold(hca,'off');
    hca.YLabel.String = 'f [Hz]';
    irf_legend(hca,{'f_{filter}','f_{LH}'},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end
if 1 % frequencies 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,mva_l123); 
    hca.YLabel.String = 'l_{MVA,E}';
    irf_legend(hca,{'max','inter','min'},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end
if 1 % velocitiess 
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,vs);     
    clim = 0.6; 
    ind_plot = find(cmaxs(:,2)>clim);
    hold(hca,'on')
    irf_plot(hca,vs(ind_plot,:),'r*');
    hca.YLabel.String = 'v [km/s]';
    irf_legend(hca,{'','',['v(C>' num2str(clim) ')']},[0.98, 0.95]);    
    irf_zoom(hca,'y');
end

title(h(1),['nTlh = ' num2str(nTlh) ', f_{filt} = ' num2str(ffiltlim) 'f_{LH}'])
irf_zoom(h,'x',tintlim)
irf_plot_axis_align

%% Get smaller scale images
ind_tints = [40 41 42 43 44 69 73 74 75 93 94 95 97];
for itint = ind_tints
    tint = tints(itint,2:3)+0.02*[-1 1];
    flim = 0.5;
    v=linspace(20,800,30);
    
    % Obtain correlation parameters
    tool.single_event
    
    irf_match_phibe_vis('best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir);
    pause
    close
end

tintlim = toepoch([2007 04 17 15 32 19.86;2007 04 17 15 32 19.82])'; % C3#1

%% do direction gif on longer time interval
tint = toepoch([2007 04 17 15 32 19;2007 04 17 15 32 24])'; % C3#1
%%
sc = 3;
flim = 0.5;
v=linspace(20,800,30);

% Obtain correlation parameters
tool.single_event

irf_match_phibe_vis('best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir);

%% Visualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass);      
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_dir_.gif'],'DelayTime',0.01,'LoopCount',inf);
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);

%% Add the event and save data to TimeTable
comment = 'Small fluxrope.'
toSave = 'single';
tool.add_event

%%
tint = irf_time(['2007-04-17T15:32:19.411252975Z';'2007-04-17T15:32:19.656185388Z'],'utc>epoch')';
tint = irf_time(['2007-04-17T15:32:19.695947170Z';'2007-04-17T15:32:19.994955539Z'],'utc>epoch')';
tint = irf_time(['2007-04-17T15:32:19.966327190Z';'2007-04-17T15:32:20.244659423Z'],'utc>epoch')';