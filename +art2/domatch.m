% Matching plot
% tint defined before
cd /Users/Cecilia/Data/BM/20070831

if 1 % Constants
    mu0=4*pi*1e-7;
    e=1.6e-19;
end

c_eval('peaTepar?z=irf_resamp(irf_tlim(parTe?,[tint(1)-15 tint(2)+15]),EparAC?);',3:4);
c_eval('Tm?=cn.mean(irf_tlim(peaTepar?z,tint));',3:4);
Teparav=mean([Tm3(:,2);Tm4(:,2)])*8.61734*10;
c_eval('peaNe?z=irf_resamp(irf_tlim(peaNe?,[tint(1)-15 tint(2)+15]),Epar?);',3:4);
c_eval('nm?=cn.mean(irf_tlim(peaNe?z,tint));',3:4);
Ne = mean([nm3(:,2) nm4(:,2)]);
%Ne = 0.2; Teparav = 3500;
ld = irf_plasma_calc(1,Ne,0,Teparav,Teparav,'Ld');
%    c_eval('peaTeper?z=cn_toepoch(t1,t2,perTe?);',3:4);
%    Teperav=mean([peaTeper3z(:,2);peaTeper4z(:,2)])*8.61734*10;
%    Teav=(Teparav+Teperav)/2;

if 1 % Do time matching, get v
    E4 = irf_tlim(EparAC4,tint);
    E3 = irf_resamp(irf_tlim(EparAC3,tint),E4);
    
    %%%%%% if its the short one
%    E4 = irf_tlim(diE4(:,[1 3]),tint);
%    E3 = irf_resamp(irf_tlim(diE3(:,[1 3]),tint),E4);
%    Teparav = 3500;
    %%%%%%
    
    window = 40;err=0.05;
    %[dtc dt_error tplus tminus corr]=cn_mycorr(facEac3(:,[1 4]),facEac4(:,[1 4]),window,450,err);
    [dtx corrx] = cn_xcorr(E3,E4,window,'x');
    dt = -dtx;
    dz = cn.mean(irf_tlim(dpar,tint),1);
    v=-dz/dt;
    %vav=sum(v)/size(v,1);
end
if 1 % Integrate potential and do shift
    c_eval('Phi?=irf_integrate(E?);',3:4)
    c_eval('Phi?(:,2)=Phi?(:,2)*v;',3:4)
    Phishift4=Phi4; Phishift4(:,1)=Phishift4(:,1)-dt;
    Eshift4=E4; Eshift4(:,1)=Eshift4(:,1)-dt;
end

c_eval('Ep? = irf_tlim(Eper?,tint);',3:4)
%c_eval('Ep? = irf_tlim(diE?(:,[1 2]),tint);',3:4)
Epershift4=Ep4; Epershift4(:,1)=Epershift4(:,1)-dt;
% Make figure
set(0,'defaultLineLineWidth', 1);
fig = figure(22);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(fig,'position',[10 416 465 463]);
set(fig,'position',[63 234 330 463]);
nPanels = 3;
h = irf_plot(nPanels);
set(h(1),'position',[0.1900    0.5600    0.700    0.23])
set(h(2),'position',[0.1900    0.3300    0.700    0.23])
set(h(3),'position',[0.1900    0.1000    0.700    0.23])

isub = 1;

if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,E3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,E4,'color','b')
    irf_plot(hca,Eshift4,'b--')
    ylabel(hca,'E_{||} [mV/m]')
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
if 1 % Perpendicular electric field
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Eper3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Eper4,'color','b')
    irf_plot(hca,Epershift4,'b--')
    ylabel(hca,'E_{\perp} [mV/m]')
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
if 1 % Parallel electrostatic potential
    hca = h(isub); isub = isub + 1;
    %irf_plot(hca,{E3,E4,E4},'comp','dt',[0 0 dt])
    irf_plot(hca,Phi3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Phi4,'color','b')
    irf_plot(hca,Phishift4,'b--')
    ylabel(hca,'\phi')
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end


irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');

% Labeling
abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.02,0.92],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:(isub-1)%nPanels
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})
    grid(h(k),'off')
    set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end

% Add a length axis on top
if 1
    ax1=h(1);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    % put each 2rd tick
    tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    tticks = xticks0(ticks);
    lticks = round(tticks*v-xticks0(tick0)*v);
    
    %nticksl=xlim(1)/5/r_e;
    %nticksr=xlim(2)/5/r_e;
    %xticks=0:xlim(2)/nticksr:xlim(2);
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end
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
if exist('ind','var') && ind==2
tticklabels = get(h(3),'xticklabel');

tticklabels{6} = ' ';
tticklabels{16} = ' ';
set(h(3),'xticklabel',tticklabels)
end

% Finishing touches to figure
titlestr = {['dt = ' num2str(abs(dt*1000),'%.0f') ' ms , v = ' num2str(v,'%.0f') ' km/s,'],...
            ['T_{e} = ' num2str(round(Teparav/100)*100,'%.0f') ' eV , \lambda_{De} = ' num2str(ld/1000,'%.1f') ' km']};
ht = title(h(1),titlestr);
%if ind==15
%    set(ht,'position',[0.16 34 17])
%elseif ind==2
    %set(ht,'position',[-6 34 17])
    %set(ht,'position',[0.16 34 17])
%end

% Fix ylims for some axes
irf_zoom(h,'y')
set(fig,'position',[63 234 330 463]);

