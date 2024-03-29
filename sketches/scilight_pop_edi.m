syms x y lx ly phi0

R = [x y];
phi = phi0*exp(-x^2/lx^2 -y^2/ly^2);

E = -gradient(phi,R);
n = divergence(E,R);

mf_phi = matlabFunction(phi,'Vars',[x,y,lx,ly,phi0]);
mf_Ex = matlabFunction(E(1),'Vars',[x,y,lx,ly,phi0]);
mf_Ey = matlabFunction(E(2),'Vars',[x,y,lx,ly,phi0]);
mf_n = matlabFunction(n,'Vars',[x,y,lx,ly,phi0]);

xv = linspace(-3,3,100);
yv = linspace(-4,4,100);
[XV,YV,ZV] = meshgrid(xv,yv,0);
y0 = 0;
lx = 1;
ly = 1.5;
phi0 = 1.1;



[XE,YE] = meshgrid(linspace(-3,3,21),linspace(-3,3,21));
EX = mf_Ex(XE,YE,lx,ly,phi0);
EY = mf_Ey(XE,YE,lx,ly,phi0);

cmap = pic_colors('blue_white');
blue = cmap(fix(size(cmap,1)/2),:);

hca = subplot(2,1,1); h(1) = hca;
hl = plot(hca,xv,-mf_phi(xv,yv,lx,ly,phi0),xv,mf_Ex(xv,yv,lx,ly,phi0),'k',xv,-mf_n(xv,yv,lx,ly,phi0));
hl(1).Color = blue;
hl(2).Color = [0 0 0];
c_eval('hca.Children(?).LineWidth = 1;',1:3)

% hca = subplot(3,1,2);
% pcolor(hca,yv,xv,mf_phi(XV,YV,lx,ly,phi0))
% colormap(hca,flipdim(pic_colors('blue_white'),1))
% shading(hca,'flat')
% hold(hca,'on')
% quiver(hca,XE,YE,EX,EY,'k')
% hold(hca,'off')

hca = subplot(2,1,2); h(2) = hca;
pcolor(hca,xv,yv,mf_phi(XV,YV,lx,ly,phi0)')
colormap(hca,flipdim(pic_colors('blue_white'),1))
shading(hca,'flat')
hold(hca,'on')
quiver(hca,XE,YE,EX,EY,'k')
hold(hca,'off')
%axis(hca,'equal')

hlinks = linkprop(h(1:2),{'XLim'});

%%
hca = subplot(3,1,3);
%surf(hca,XV,YV,mf_n(XV,YV,lx,ly,phi0),mf_phi(XV,YV,lx,ly,phi0))
surf(hca,XV,YV,mf_phi(XV,YV,lx,ly,phi0),mf_n(XV,YV,lx,ly,phi0))
colormap(hca,pic_colors('blue_red'))
hca.CLim = [-3 3];
colorbar('peer',hca)
axis(hca,'equal');
shading(hca,'flat');
%%

[XE,YE] = meshgrid(linspace(-3,3,15),linspace(-3,3,17));
EX = mf_Ex(XE,YE,lx,ly,phi0);
EY = mf_Ey(XE,YE,lx,ly,phi0);

cmap = pic_colors('blue_white');
blue = cmap(fix(size(cmap,1)/2),:);
colors = pic_colors('matlab');

hca = subplot(3,1,[1 2]); h(1)d = hca;
pcolor(hca,xv,yv,mf_phi(XV,YV,lx,ly,phi0))
colormap(hca,flipdim(pic_colors('blue_red'),1))
shading(hca,'flat')
hold(hca,'on')
quiver(hca,XE,YE,EX,EY,'k')
hold(hca,'off')
%hca.Position(2) = 0.4;
%hca.Position(4) = 0.5;

hca = subplot(3,1,3); h(2) = hca;
hl = plot(hca,xv,mf_phi(xv,0,lx,ly,phi0),xv,mf_Ex(xv,0,lx,ly,phi0),'k',xv,-0.5*mf_n(xv,0,lx,ly,phi0));
hl(3).Color = blue;
hl(2).Color = [0 0 0];
hl(1).Color = colors(2,:);
c_eval('hl(?).LineWidth = 4;',1:3)
%c_eval('axis(h(?),''off'');',1:2)

%hlinks = linkprop(h(1:2),{'XLim'});

%%
%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

twpe = 24000;
xlim = [80 140];
zlim = [-7 0];

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

pic.plot_map({'tepar'},'cmap',{pic_colors('candy4')},'clim',{[0 0.3]},'sep','smooth',1);
%pic.plot_map({'Epar'},'cmap',{pic_colors('blue_red')},'clim',{[-1 1]},'sep','smooth',2);

%% First run paper_bgk_timeseries

%%
%fig = figure;
npanels = 2;
h = irf_plot(npanels);
clear h_all;
h_all = h;
isub = 1;
fclim = [0 0.0025];

vlim_f = [-29 15];
ylim1 = [0 3.6];
ylim2 = [0 360];
  
   
  if 1 % F
    %hca = h(isub); isub = isub + 1;
    hca = irf_panel('fmodel map');
    i_vlim = find((Fspecrec.f>(vlim_f(1)*1.01)).*(Fspecrec.f<(vlim_f(2)*1.01)));
    Fspecrec_plot = Fspecrec;
    Fspecrec_plot.p = Fspecrec_plot.p(:,i_vlim);
    Fspecrec_plot.f = Fspecrec_plot.f(i_vlim);
    irf_spectrogram(hca,Fspecrec_plot,'lin');
    hca.YLabel.String = {'v_{||}','(10^3 km/s)'};
    edi_color = [1 1 1];
    vph_color = [1 1 1];
    h_all_markings = [];
    if 1 % EDI energies
      hold(hca,'on')
      %lines_EDI_plus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      lines_EDI_minus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      %hleg_EDI = irf_legend(hca,{'- EDI'},[0.08 0.99],'color',lines_EDI_minus(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_EDI_plus; lines_EDI_minus; hleg_EDI];
      %h_all_markings = [h_all_markings; lines_EDI_minus; hleg_EDI]; 
    end
    if 0 % model phase velocity
      hold(hca,'on')
      line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
      lines_vphav = irf_plot(hca,tsVph*1e-6,'LineWidth',1.0,'Color',vph_color,'LineStyle','--');
      %lines_vphav = irf_plot(hca,tsVph*1e-6,'--k');
      %hleg_vphav = irf_legend(hca,{'-- v_{ph,mod}'},[0.32 0.99],'color',lines_vphav(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_vphav; hleg_vphav];
    end
    if 0 % observed phase velocity
      hold(hca,'on')
      lines_vphobs = irf_plot(hca,tsVphIndividual*1e-3,'*k','LineWidth',1.5,'Color',vph_color);
      lines_vphobs.MarkerSize = 4;
      %hleg_vphobs = irf_legend(hca,{'* v_{ph,obs}'},[0.18 0.99],'color',lines_vphobs(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_vphobs; hleg_vphobs];
    end  
    colormap(hca,cn.cmap('white_blue'))
    if 0 % str info
      str_info = {'unperturbed f:';...
        ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
        ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
        ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
        sprintf('beta_{Schamel}=%g',beta);...
        };
      irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
    end
    hca.YLim = vlim_f;
    %hca.CLim = fclim;
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.Interpreter = 'tex';
  end
 
  if 1 % edi, phi comparison for one spacecraft
    isub = isub + 1;
    hca = irf_panel('edi phi comp');  
    set(hca,'ColorOrder',colors)
    palim = 180;
    c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6},''color'',colors(1,:));',ic)
    hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};  
    hca.YLabel.Interpreter = 'tex';

    set(hca,'ColorOrder',colors)

    hca.YLim = ylim1; 
    yyaxis('right')
    ax = gca;
    irf_plot(ax,phi1)
    ax.YLabel.String = {'\phi (V)'};
    ax.YLabel.Interpreter = 'tex';
    ax.YLim = ylim2;
    %ax.YLabel.Position(1) = 1.07;
    %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.02 0.99],'fontsize',12,'color',[0 0 0]);
    %hca.Title.String = 'mms 1';
    %hca.YLabel.String = 'Electron flux';
    %ax.YLabel.String = 'Electrostatic potential';
    hca.YAxis(1).Label.String = {'Electron flux','(10^6 s^{-1}cm^{-2}sr^{-1})'};
    hca.YAxis(2).Label.String = {'Electrostatic potential','(V)'};
  end

irf_zoom(h,'x',tint_figure)
  irf_plot_axis_align([h ax])
  %%
  h(3).YLabel.String = {'v_{||}','(10^3 km/s)'};
  h(3).YLabel.Interpreter = 'tex';
  
  h(2).YLabel.String = {'\phi','(V)'};  
  h(2).YLabel.Interpreter = 'tex';
  h(2).YLabel.Position = h(1).YLabel.Position;
  
  iref = 1;
  ijmp = 1;
  %h_all(iref+4+ijmp).Position = h_all(iref+4).Position;  
  
  
  %h_all(iref+7).Position = h_all(iref+5).Position;
  
  %irf_zoom(h_all,'x',[tint_phi(1) ts_edi_flux180.time(end)])  % tint_phi
  %irf_zoom(h_all,'x',[tint_model(1) tint_model(end)])  % tint_phi
  irf_zoom(h_all,'x',tint_figure)  % tint_phi
  
  h_all(iref+4+ijmp).XLabel = [];
  
  h(1).YLim = [-70 70];
  h(4).YLim = [-0.9 1.7];
  
  ax2_flux.Position = ax1_flux.Position;
  ax2_dn.Position = ax1_dn.Position;
  ax2_dn.YTick = ax1_dn.YTick/1e-3/n0;
  ax2_dn.YLim = ax1_dn.YLim/1e-3/n0;

  
  h(3).YLabel.Position(1) = h(1).YLabel.Position(1);
  hold(h(1),'on')
  plot(h(1),h(1).XLim,[0 0],'color',[0.8 0.8 0.8])
  hold(h(1),'off')
  hold(h(3),'on')
  plot(h(3),h(3).XLim,[0 0],'color',[0.6 0.6 0.6])
  hold(h(3),'off')

  %h_all(iref+7).XLabel = [];
  %h_all(end-2).YAxisLocation = 'right';
  %h_all(end-3).YAxisLocation = 'right';
  legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
  legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
  for ipanel = 1:npanels
    irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  end
  c_eval('h_all(?).XGrid = ''off''; h_all(?).YGrid = ''off'';',1:numel(h_all))
  
  h(3).YLim = [-27 15];
  
  % annotation
  h(3).YLim = vlim_f;
  irf_zoom(h,'x',tint_figure)
  h(3).CLim = [0 2.5]*1e-3;
  h(3).YLim = [-27 15];
  
  legend(h(3).Children([3 4]),{'v_{ph}','v_{ph,av}'},'box','off','Orientation','Horizontal')
  annotation('textarrow',[0.685 0.685]-0.025,[0.47 0.495]+0.00,'string',{'EDI range'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.76 0.76]-0.015,[0.550 0.500]+0.02,'string',{'v_{ph,av}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.30 0.325]+0.015,[0.550 0.510]+0.02,'string',{'v_{ph,ind}'},'fontsize',12,'horizontalalignment','center');
  