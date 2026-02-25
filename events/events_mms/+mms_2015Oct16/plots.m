%% Figure: Parallel electric field and parallel electron velocity and energy dissipation
tint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere

h = irf_plot(15);
ic = 4;
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % ePDist omni 32
  hca = irf_panel('e DEF omni 32');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist omni 32
  hca = irf_panel('e DEF omni 64');
  c_eval('irf_spectrogram(hca,ePDist?.e64.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.pitchangles(dmpaB?,[]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 0 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).x,mvaE2.tlim(tint).x,mvaE3.tlim(tint).x,mvaE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).y,mvaE2.tlim(tint).y,mvaE3.tlim(tint).y,mvaE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EL par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1par.tlim(tint),mvaE2par.tlim(tint),mvaE3par.tlim(tint),mvaE4par.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EL perp
  hca = irf_panel('EL perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).x,mvaE2perp.tlim(tint).x,mvaE3perp.tlim(tint).x,mvaE4perp.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_{\perp,L}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EM perp
  hca = irf_panel('EM perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).y,mvaE2perp.tlim(tint).y,mvaE3perp.tlim(tint).y,mvaE4perp.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_{\perp,M}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN perp
  hca = irf_panel('EN perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1perp.tlim(tint).z,mvaE2perp.tlim(tint).z,mvaE3perp.tlim(tint).z,mvaE4perp.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_{\perp,N}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0
%load('caa/cmap.mat');
for ii = [4 5 7 8]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:10
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
end
end
irf_zoom(h,'x',tint)
%irf_zoom(h([1:3 6 10]),'y')
irf_plot_axis_align

%% Figure: Plot all length scales and frequencies

h = irf_plot(14);
if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % vtp
  hca = irf_panel('vt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{vtp1.tlim(tint),vtp2.tlim(tint),vtp3.tlim(tint),vtp4.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{tp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % vte
  hca = irf_panel('vte');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{vte1.tlim(tint),vte2.tlim(tint),vte3.tlim(tint),vte4.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{te}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % vA
  hca = irf_panel('vA');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{vA1.tlim(tint),vA2.tlim(tint),vA3.tlim(tint),vA4.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{A}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % flh
  hca = irf_panel('flh');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{flh1.tlim(tint),flh2.tlim(tint),flh3.tlim(tint),flh4.tlim(tint)},'comp');
  hca.YLabel.String = {'f_{LH}','(Hz)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % fce
  hca = irf_panel('fce');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{fce1.tlim(tint),fce2.tlim(tint),fce3.tlim(tint),fce4.tlim(tint)},'comp');
  hca.YLabel.String = {'f_{ce}','(Hz)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % fpe
  hca = irf_panel('fpe');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{fpe1.tlim(tint),fpe2.tlim(tint),fpe3.tlim(tint),fpe4.tlim(tint)},'comp');
  hca.YLabel.String = {'f_{pe}','(Hz)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % fpp
  hca = irf_panel('fpp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{fpp1.tlim(tint),fpp2.tlim(tint),fpp3.tlim(tint),fpp4.tlim(tint)},'comp');
  hca.YLabel.String = {'f_{pi}','(Hz)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end

if 1 % Le
  hca = irf_panel('Le');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Le1.tlim(tint),Le2.tlim(tint),Le3.tlim(tint),Le4.tlim(tint)},'comp');
  hca.YLabel.String = {'L_{e}','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
end
if 1 % Ld
  hca = irf_panel('Ld');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Ld1.tlim(tint),Ld2.tlim(tint),Ld3.tlim(tint),Ld4.tlim(tint)},'comp');
  hca.YLabel.String = {'L_{d}','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
end
if 1 % re
  hca = irf_panel('re');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{re1.tlim(tint),re2.tlim(tint),re3.tlim(tint),re4.tlim(tint)},'comp');
  hca.YLabel.String = {'r_{e}','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
end
if 1 % rp
  hca = irf_panel('rp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{rp1.tlim(tint),rp2.tlim(tint),rp3.tlim(tint),rp4.tlim(tint)},'comp');
  hca.YLabel.String = {'r_{p}','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
end
if 1 % Lp
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Lp1.tlim(tint),Lp2.tlim(tint),Lp3.tlim(tint),Lp4.tlim(tint)},'comp');
  hca.YLabel.String = {'L_{p}','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
end  