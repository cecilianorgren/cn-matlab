%% Flux transport

tref = EpochTT('2017-07-25T22:09:00.00Z'); % Maximum |B|
c_eval('intExperp? = irf_integrate(lmnE?perp.x,tref);')
c_eval('intEyperp? = irf_integrate(lmnE?perp.y,tref);')
c_eval('intEzperp? = irf_integrate(lmnE?perp.z,tref);')

dAy = 23*1e-3; % from intEyperp between minimum (which is close to max B) and just before DF


%% Observations
ic = 1;
npanels = 6;
h = irf_plot(npanels);
cmap = colormap(pic_colors('candy4'));
fontsize = 12;

t1 = EpochTT('2017-07-25T22:09:45.3Z');
t2 = EpochTT('2017-07-25T22:10:06.0Z');

tint_figure = irf.tint('2017-07-25T22:09:00.000Z/2017-07-25T22:10:18.000Z');

if 1 % B LMN 
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?.x,lmnB?.y,lmnB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end 
if 1 % E LMN
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  fhigh = 2;
  c_eval('irf_plot(hca,{lmnE?perp.x.filt(0,fhigh,[],3),lmnE?perp.y.filt(0,fhigh,[],3),lmnE?perp.z.filt(0,fhigh,[],3)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp} (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_legend(hca,{sprintf('f<%g Hz',fhigh)},[0.05 0.98],'fontsize',fontsize,'color',[0 0 0]);
end 
if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?perp.x,lmnVi?perp.y,lmnVi?perp.z},''comp'');',ic)    
  hca.YLabel.String = {'v_{i\perp} (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 1 % Ey 1234
  hca = irf_panel('Ey');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{lmnE1perp.y,lmnE2perp.y,lmnE3perp.y,lmnE4perp.y},'comp');
  hca.YLabel.String = {'E_{y\perp} (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
end
if 1 % Ey 1234
  hca = irf_panel('Ey filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{lmnE1perp.y.filt(0,fhigh,[],3),lmnE2perp.y.filt(0,fhigh,[],3),lmnE3perp.y.filt(0,fhigh,[],3),lmnE4perp.y.filt(0,fhigh,[],3)},'comp');
  hca.YLabel.String = {'E_{y\perp} (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{sprintf('f<%g Hz',fhigh)},[0.05 0.98],'fontsize',fontsize,'color',[0 0 0]);
end
if 1 % int Ey dt
  hca = irf_panel('int Ey dt LMN');
  fhigh = 2;
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{intEyperp1,intEyperp2,intEyperp3,intEyperp4},'comp');
  hca.YLabel.String = {'\int E_{y\perp} dt (mV/m s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
end

irf_zoom(h,'x',tint_figure)
irf_zoom(h,'y')

irf_pl_mark(h,t1,'black','linestyle','--')
irf_pl_mark(h,t2,'black','linestyle','--')

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.0;',1:numel(hl))

c_eval('h(?).FontSize = fontsize;',1:numel(h))


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:numel(h)
  irf_legend(h(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;  
end



%% Harris current sheet
units = irf_units;

Bx = @(z,B0,L) B0*tanh(z/L);
Ay = @(z,B0,L) B0*L*log(cosh(z/L));
n = @(z,n0,L) n0*cosh(z/L).^(-2);

n0 = 0.3*1e6; % mm
B0 = 20*1e-9; % T
L = 2000*1e3; % m
zvec = 2*L*linspace(-1,1,150);
i_pre_df = find(Ay(zvec,B0,L)<dAy);
vA = Bx(zvec,B0,L)./sqrt(units.mu0*units.mp*n(zvec,n0,L));

fontsize = 12;

h = setup_subplots(4,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zvec/L,Bx(zvec,B0,L))
hca.XLabel.String = 'z/L';
hca.YLabel.String = 'B_x (T)';

hca = h(isub); isub = isub + 1;
plot(hca,zvec/L,Ay(zvec,B0,L))
hold(hca,'on')
plot(hca,zvec/L,zvec*0+dAy)
hold(hca,'off')
hca.XLabel.String = 'z/L';
hca.YLabel.String = 'A_y (Tm)';

hca = h(isub); isub = isub + 1;
plot(hca,zvec/L,n(zvec,n0,L)*1e-6)
hca.XLabel.String = 'z/L';
hca.YLabel.String = ['n (cc)'];

hca = h(isub); isub = isub + 1;
plot(hca,zvec/L,abs(vA)*1e-3)
hca.XLabel.String = 'z/L';
hca.YLabel.String = 'v_A (km/s)';


irf_legend(h(1),sprintf('L = %.0f km',L*1e-3),[0.02 0.98],'fontsize',fontsize,'color','k')

for ip = 1:numel(h)
  hold(h(ip),'on')
  plot(h(ip),zvec([i_pre_df([1 1])])/L,h(ip).YLim,'k--')
  plot(h(ip),zvec([i_pre_df([end end])])/L,h(ip).YLim,'k--')

  %plot(h(ip),[1 1]*L*1e-3,h(ip).YLim,'color',[0    0.4470    0.7410])
  %plot(h(ip),-1*[1 1]*L*1e-3,h(ip).YLim,'color',[0    0.4470    0.7410])
  
  hold(h(ip),'off')
end


hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.0;',1:numel(hl))

c_eval('h(?).FontSize = fontsize;',1:numel(h))
