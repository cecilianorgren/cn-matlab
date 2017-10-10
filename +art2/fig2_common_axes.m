% Matching plot

cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields
% 5: nice mono also, 2: ugly mono
tints{1} = tint{highQuality(2)};
tints{2} = tint{highQuality(5)};

mu0=4*pi*1e-7;
e=1.6e-19;
getData = 0;
if getData
    c_eval('peaTepar?z=irf_resamp(irf_tlim(parTe?,[tint(1)-15 tint(2)+15]),EparAC?);',3:4);
    c_eval('Tm?=cn.mean(irf_tlim(peaTepar?z,tint));',3:4);
    Teparav=mean([Tm3(:,2);Tm4(:,2)])*8.61734*10;
    c_eval('peaNe?z=irf_resamp(irf_tlim(peaNe?,[tint(1)-15 tint(2)+15]),Epar?);',3:4);
    c_eval('nm?=cn.mean(irf_tlim(peaNe?z,tint));',3:4);
    Ne = mean([nm3(:,2) nm4(:,2)]);
end
% Average of two concerned intervals
Ne = 0.04; Teparav = 1600;
ld = irf_plasma_calc(1,Ne,0,Teparav,Teparav,'Ld');

%
if 1 % Do time matching, get v
    c_eval('E4? = irf_tlim(Epar4,tints{?});',1:2);
    c_eval('E3? = irf_resamp(irf_tlim(Epar3,tints{?}),E4?);',1:2);
    
    %%%%%% if its the short one
%    E4 = irf_tlim(diE4(:,[1 3]),tint);
%    E3 = irf_resamp(irf_tlim(diE3(:,[1 3]),tint),E4);
%    Teparav = 3500;
    %%%%%%
    
    window = 40;err=0.05;
    %[dtc dt_error tplus tminus corr]=cn_mycorr(facEac3(:,[1 4]),facEac4(:,[1 4]),window,450,err);
    c_eval('[dtx? corrx?] = cn_xcorr(E3?,E4?,window,''x'');',1:2);
    c_eval('dt? = -dtx?;',1:2);
    c_eval('dz? = cn.mean(irf_tlim(dpar,tints{?}),1);',1:2);
    c_eval('v?=-dz?/dt?;',1:2);
    %vav=sum(v)/size(v,1);
    
if 1 % Integrate potential and do shift
    c_eval('Phi?!=irf_integrate(E?!);',3:4,1:2);
    c_eval('Phi?!(:,2)=Phi?!(:,2)*v!;',3:4,1:2);
    c_eval('Phishift?!=Phi?!; Phishift?!(:,1)=Phishift?!(:,1)-dt!;,',3:4,1:2);
    c_eval('Eshift?!=E?!; Eshift?!(:,1)=Eshift?!(:,1)-dt!;',3:4,1:2);
end
end

c_eval('Ep?! = irf_tlim(Eper?,tints{!});',3:4,1:2)
%c_eval('Ep? = irf_tlim(diE?(:,[1 2]),tint);',3:4)
c_eval('Epershift?!=Ep?!; Epershift?!(:,1)=Epershift?!(:,1)-dt!;',3:4,1:2);
%
% Make figure
set(0,'defaultLineLineWidth', 1);
fig = figure(22);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(fig,'position',[10 416 465 463]);
set(fig,'position',[563 234 890 700]);
nPanels = 6;
h = irf_plot(nPanels);
ww = 0.3; hh = 0.23; x1 = 0.19; x2 = 0.50;
set(h(1),'position',[x1 0.5600    ww hh])
set(h(2),'position',[x1 0.3300    ww hh])
set(h(3),'position',[x1 0.1000    ww hh])
set(h(4),'position',[x2 0.5600    ww hh])
set(h(5),'position',[x2 0.3300    ww hh])
set(h(6),'position',[x2 0.1000    ww hh])
%
isub = 1;

if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,E31,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,E41,'color','b')
    irf_plot(hca,Eshift41,'b--')
    ylabel(hca,'E_{||} [mV/m]')
    %set(hca,'colororder',[0.2 0.5 0.2;0 0 1])
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
end
%
if 1 % Perpendicular electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Eper3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Eper4,'color','b')
    irf_plot(hca,Epershift41,'b--')
    ylabel(hca,'E_{\perp} [mV/m]')
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
end
%
if 1 % Parallel electrostatic potential
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Phi31,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Phi41,'color','b')
    irf_plot(hca,Phishift41,'b--')
    ylabel(hca,'\phi_{||} [V]')
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
    if 0 % Ytick and Yticklabel ePhi/Te
        %%
        %set(gca,'Box','off');
        axesPosition=get(hca,'Position');
        yticks0=get(hca,'ytick');
        yticklabels_num=yticks0'/Teparav;
        yticklabels_str=num2str(yticklabels_num,'%.2f');
        ylim0=get(hca,'ylim');
        xlim0=get(hca,'xlim');
        ylim=ylim0/Teparav;
        xticks0=get(hca,'xtick');
        ax2=axes('Position',axesPosition,'ylimmode','manual','yaxislocation','right',...
            'color','none','xaxislocation','top','xtick',[],...
            'xticklabel',[],'ytick',yticks0,'yticklabel',yticklabels_str,...
            'Box','off','GridLineStyle','none',...
            'ylimmode','manual','ylim',ylim0,'xlim',xlim0);
        ylabel(ax2,'e\phi/ T_{e}')
    end
end
%
if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,E32,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,E42,'color','b')
    irf_plot(hca,Eshift42,'b--')
    ylabel(hca,'E_{||} [mV/m]')
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
end
if 1 % Perpendicular electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Eper3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Eper4,'color','b')
    irf_plot(hca,Epershift42,'b--')
    ylabel(hca,'E_{\perp} [mV/m]')
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
end
if 1 % Parallel electrostatic potential
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Phi32,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Phi42,'color','b')
    irf_plot(hca,Phishift42,'b--')
    ylabel(hca,'\phi [V]')
    %irf_legend(hca,{'C3','C4'},[0.98 0.9])
    irf_zoom(hca,'y')
    
end
%

irf_zoom(h(1:3),'x',tints{1});
irf_zoom(h(4:6),'x',tints{2});
set(h([1 4]),'ylim',[-25 25])
set(h([2 5]),'ylim',[-4 19])
set(h([3 6]),'ylim',[-40 390])%set(h([3 6]),'ylim',[-70 240])
%irf_plot_axis_align;
irf_timeaxis(hca,'usefig');
for ff = 4:6
    set(h(ff),'yticklabel',[])
    ylabel(h(ff),' ')
end
tlabels1 = get(h(3),'xticklabel'); tlabels1{1} = ' '; tlabels1{11} = ' '; tlabels1{end} = ' ';
tlabels2 = get(h(6),'xticklabel'); tlabels2{1} = ' ';
set(h(6),'xticklabel',tlabels2)
set(h(3),'xticklabel',tlabels1)

if 1 % Ytick and Yticklabel ePhi/Te
        %%
        %set(gca,'Box','off');
        hca = h(6);
        axesPosition=get(hca,'Position');
        yticks0=get(hca,'ytick');
        yticklabels_num=yticks0'/Teparav;
        yticklabels_str=num2str(yticklabels_num,'%.2f');
        ylim0=get(hca,'ylim');
        xlim0=get(hca,'xlim');
        ylim=ylim0/Teparav;
        xticks0=get(hca,'xtick');
        ax2=axes('Position',axesPosition,'ylimmode','manual','yaxislocation','right',...
            'color','none','xaxislocation','top','xtick',[],...
            'xticklabel',[],'ytick',yticks0,'yticklabel',yticklabels_str,...
            'Box','off','GridLineStyle','none',...
            'ylimmode','manual','ylim',ylim0,'xlim',xlim0);
        ylabel(ax2,'e\phi/ T_{e}')
    end
%

% Labeling
abcde={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.04,0.92],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:3%nPanels
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})    
    %set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    %irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end
for k=4:6%nPanels
    %irf_legend(h(k),x[abcde(k)],legloc{1},'color',colors{1})
    
    set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end
for k=1:(isub-1)%nPanels
    %irf_legend(h(k),x[abcde(k)],legloc{1},'color',colors{1})
    grid(h(k),'off')
    %set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    %irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end
%
% Add a length axis on top
if 1
    %
    % Left panels
    %
    ax1=h(1);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    % put each 2rd tick
    tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    tticks = xticks0(ticks);
    lticks = round(tticks*v1-xticks0(5)*v1);
    %lticks = round(tticks*v1-xticks0(tick0)*v1);
    
    %nticksl=xlim(1)/5/r_e;
    %nticksr=xlim(2)/5/r_e;
    %xticks=0:xlim(2)/nticksr:xlim(2);
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end    
    xticklabels{1} = ' '; 
    xticklabels{2} = ' ';
    xticklabels{end} = ' '; 
    xticklabels{end-1} = 'km'; 
    yticks=[];
    yticklabels=[];
    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
             'XAxisLocation','top','YAxisLocation','right',...
             'Color','none','Ytick',yticks,...
             'YTickLabel',yticklabels,...
             'xlim',xlim,'xtick',tticks,...
             'XTickLabel',xticklabels);
    
    %
    % Right panels
    %
    ax1=h(4);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    % put each 2rd tick
    tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    ticks = ticks(1:end-1) + 1;
    tticks = xticks0(ticks);
    lticks = round(tticks*v2-xticks0(3)*v2);
    %lticks = round(tticks*v2-xticks0(tick0)*v2);
    
    %nticksl=xlim(1)/5/r_e;
    %nticksr=xlim(2)/5/r_e;
    %xticks=0:xlim(2)/nticksr:xlim(2);
    xticklabels = cell(1,numel(lticks));
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end
    %xticklabels{1} = ' '; 
    %xticklabels{end} = ' '; 
    xticklabels{end} = 'km'; 
    
    yticks=[];
    yticklabels=[];
    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
             'XAxisLocation','top','YAxisLocation','right',...
             'Color','none','Ytick',yticks,...
             'YTickLabel',yticklabels,...
             'xlim',xlim,'xtick',tticks,...
             'XTickLabel',xticklabels);         
end
%
%if exist('ind','var') && ind==2
%tticklabels = get(h(3),'xticklabel');

%tticklabels{6} = ' ';
%tticklabels{16} = ' ';
%set(h(3),'xticklabel',tticklabels)
%end
%%

% Finishing touches to figure
titlestr = {['dt = ' num2str(abs(dt*1000),'%.0f') ' ms , v = ' num2str(v,'%.0f') ' km/s,'],...
            ['T_{e} = ' num2str(round(Teparav/100)*100,'%.0f') ' eV , \lambda_{De} = ' num2str(ld/1000,'%.1f') ' km'],...
            [' ']};
ht = title(h(1),titlestr);
if ind==15
    set(ht,'position',[0.16 34 17])
elseif ind==2
    %set(h(3),'ylim',get(h(3),'ylim')+[-5 5])
    %set(ht,'position',[-6 34 17])
    %set(ht,'position',[0.16 34 17])
end

% Fix ylims for some axes
%irf_zoom(h(1:2),'y')
set(fig,'position',[63 234 370 463]);

%irf_plot_axis_align;

