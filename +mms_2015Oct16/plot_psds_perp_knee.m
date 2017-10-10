%%
% Automatic results made from gui above
ic = 1:4;
tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
c_eval('[out?,l?,mva?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
mva = mva3;
% rotate around x-axis 
turn_angle = 00; % degrees
tmpv = mva;
mva(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
mva(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
% v is 3x3 matrix:
%v1: [0.0906 0.0896 0.9918] - first row
%v2: [0.3510 -0.9349 0.0524] - second row
%v3: [0.9320 0.3434 -0.1161] - third row  
c_eval('mvaB?=irf.ts_vec_xyz(dmpaB?brst.time,[dmpaB?brst.dot(mva(1,:)).data_ dmpaB?brst.dot(mva(2,:)).data_ dmpaB?brst.dot(mva(3,:)).data_]);',ic)
c_eval('mvaE?=irf.ts_vec_xyz(dslE?brst.time,[dslE?brst.dot(mva(1,:)).data_ dslE?brst.dot(mva(2,:)).data_ dslE?brst.dot(mva(3,:)).data_]);',ic)
mvaJ=irf.ts_vec_xyz(jbrst.time,[jbrst.dot(mva(1,:)).data_ jbrst.dot(mva(2,:)).data_ jbrst.dot(mva(3,:)).data_]);
%% Initialize plot
[h1,h2] = mms.initialize_comined_plot(8,5,4,5,'vertical');

%tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintDist = irf.tint('2015-10-16T10:33:45.00Z',0.1);
timeOverview = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tint = timeOverview;
ic = 1;

%% Plot time series

hca = irf_panel('B_L');
irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{L}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_M');
irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{M}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('B_N');
irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{N}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Babs');
irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
hca.YLabel.String = {'|B|','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,{mvaJ.x,mvaJ.y,mvaJ.z,jbrst.dot(avB.resample(jbrst)/avB.abs)},'comp');
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J (4 SC)','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_L','J_M','J_N','J_{||}'},[0.95 0.95]);

hca = irf_panel('mvaE N');
c_eval('irf_plot(hca,mvaE?);',ic);
hca.YLabel.String = {'E','[mV/m]'};
irf_legend(hca,{'E_L','E_M','E_N'},[0.95 0.95]);
irf_legend(hca,{irf_ssub('mms?',ic)},[0.05 0.95]);

hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst.tlim(tint).tlim(tint));',ic);
hca.YLabel.String = {'v_e','(km/s)'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
hold(hca,'on');
c_eval('Tlines = irf_plot(hca,Te?_lowres.tlim(tint));',ic)
Tlines(1).Color = [0 0 0];
Tlines(2).Color = [1 0 0];
Tlines(3).Color = [0 0.7 0];
%c_eval('irf_plot(hca,veeV?)',ic)
hold(hca,'off');
irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.02 0.6],'color','cluster')
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%caxis(hca,[-15 -5])
hca.YLim = [10 30e3];
ylabel(hca,{'E_{e,OMNI}','(eV)'},'fontsize',12,'Interpreter','tex');
hca.YLim = [10 2e3];
colormap(hca,'gray')


if 0
hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {'B_{DMPA}','(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('brst E');
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E_{DSL}','(mV/m)'};
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

if 0
hca = irf_panel('brst scPot');
c_eval('irf_plot(hca,(-1)*P?brst.tlim(tint));',ic);
hca.YLabel.String = {'-scPot','(V)'};
end

hca = irf_panel('brst n');
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
c_eval('irf_plot(hca,ne?_lowres.tlim(tint));',ic);
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
irf_legend(hca,{'n_e','n_i'},[0.95 0.95]);

hca = irf_panel('brst Te');
c_eval('irf_plot(hca,Te?_lowres.tlim(tint));',ic);
hca.YLabel.String = {'T','(eV)'};
hca.YScale = 'lin';
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.95 0.95]);

hca = irf_panel('brst ve');
c_eval('irf_plot(hca,ve?brst.tlim(tint).tlim(tint));',ic);
hca.YLabel.String = {'v_e','(km/s)'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
hold(hca,'on');
c_eval('irf_plot(hca,Te?_lowres.tlim(tint))',ic)
%c_eval('irf_plot(hca,veeV?)',ic)
hold(hca,'off');
irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.5 0.1])
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%caxis(hca,[-15 -5])
hca.YLim = [10 30e3];
ylabel(hca,{'E_{e,OMNI}','(eV)'},'fontsize',12,'Interpreter','tex');

% Electron velocity in eV
me = 9.1094e-31;
e = 1.6022e-19;
c_eval('veeV? = ve?brst.abs*ve?brst.abs*1e6*0.5*me/e;',ic);
c_eval('Bhat = dmpaB?brst.resample(veeV?.tlim(tint).time); Bhat = Bhat/Bhat.abs',ic);
c_eval('veeV?par = ve?brst.dot(Bhat).times(ve?brst.dot(Bhat))*1e6*0.5*me/e;',ic)
hca=irf_panel('edistparperp');
c_eval('irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
hold(hca,'on')
c_eval('irf_plot(hca,veeV?par)',ic)
hold(hca,'off')
irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,1*[-1 1])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E','(eV)'});
colormap(hca,cn.cmap('bluered3'))

hca=irf_panel('edistparapar');
c_eval('irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,1*[-1 1])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E','(eV)'});
colormap(hca,cn.cmap('bluered3'))
end

%h1(1).Title.String = ['MMS ' num2str(ic)];
irf_zoom(h1,'x',tint)
irf_plot_axis_align
irf_zoom(h1(1:5),'y')
%for ii = 6:8
%  h1(ii).YLim = [10 500];
%end

irf_zoom(h1,'x',timeOverview)
irf_plot_axis_align
irf_zoom(h1(1:7),'y')  
hmark = 1:5;
  
  %% Check knee distributions

tint = tintDist;
vlim = 15*1e3;
elevlim = 30;
correctBin = 0;
irf_pl_mark(h1,tint.epochUnix','green')

isub = 1;
for ic = 1:4
  c_eval('dist = desDist?;',ic)

  c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
  hatVi0 = double(irf_norm(Vi0));
  c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
  hatVe0 = double(irf_norm(Ve0));

  % Get mean magnetic field direction
  c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
   

  vectors = {hatB0,'B';hatVe0,'V_e'};
  % Plot project ion onto a plane
  hca = h2(isub);  isub = isub + 1;
  vlabels = {'v_L','v_M','v_N'};
  mms.plot_projection(hca,dist,'tint',tint,'xyz',mva,'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'vlabel',vlabels,'usebincorrection',correctBin);
  hca.Title.String = {irf_ssub('MMS ?',ic),hca.Title.String{1}};
  
  hca = h2(isub); isub = isub + 1;
  vlabels = {'v_M','v_N','v_L'};
  mms.plot_projection(hca,dist,'tint',tint,'xyz',[mva(2,:); mva(3,:); mva(1,:)],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'vlabel',vlabels,'usebincorrection',correctBin);
  hca.Title.String = '';
  
  hca = h2(isub); isub = isub + 1;
  vlabels = {'v_N','v_L','v_M'};
  mms.plot_projection(hca,dist,'tint',tint,'xyz',[mva(3,:); mva(1,:); mva(2,:)],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'vlabel',vlabels,'usebincorrection',correctBin);
  hca.Title.String = '';
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; hatB0;cross(hatB0,[0 1 0])],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'usebincorrection',correctBin);
  hca.Title.String = '';
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint,'xyz',[1 0 0; cross(hatB0,[0 1 0]); hatB0],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'usebincorrection',correctBin);
  hca.Title.String = '';
end
for ii = 1:numel(h2)
  colormap(h2(ii),'jet')
  h2(ii).CLim = [-2 5];
end
h1(end).YGrid = 'off';







