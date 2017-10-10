% Matching plot

cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==1); % highest quality fields
% 5: nice mono also, 2: ugly mono
tints{1} = tint{highQuality(2)};
tints{2} = tint{highQuality(5)};


load fig_match
load Eparper
mu0=4*pi*1e-7;
e=1.6e-19;
getData = 0;
if getData
    c_eval('peaTepar?z=irf_resamp(irf_tlim(parTe?,[tint(1)-15 tint(2)+15]),Epar?);',3:4);
    c_eval('Tm?=cn.mean(irf_tlim(peaTepar?z,tint));',3:4);
    Teparav=mean([Tm3(:,2);Tm4(:,2)])*8.61734*10;
    c_eval('peaNe?z=irf_resamp(irf_tlim(peaNe?,[tint(1)-15 tint(2)+15]),Epar?);',3:4);
    c_eval('nm?=cn.mean(irf_tlim(peaNe?z,tint));',3:4);
    Ne = mean([nm3(:,2) nm4(:,2)]);
end
% Average of two concerned intervals
Ne = 0.04; Teparav = 1600;
ld = irf_plasma_calc(1,Ne,0,Teparav,Teparav,'Ld');

    
if 1 % Do time matching, get v    
    c_eval('E?! = irf_tlim(Epar?,tints{!});',3:4,1:2);    
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

% Magnetic field
% observation
%
fmin = 0.05; fmax = 0;
if 0
    c_eval('B? = irf_filt(staBpar?,fmin,fmax,5);',3:4,1:2);
    c_eval('B?! = irf_tlim(B?,tints{!});',3:4,1:2);
    %c_eval('B?! = irf_filt(B?!,fmin,fmax);',3:4,1:2);
    %h=irf_plot({B3,staBpar3})
    %irf_zoom(h,'x',tints{2})
elseif 1
    c_eval('B? = staBpar?;',3:4);
    c_eval('B?!filt = irf_filt(B?,fmin,fmax,5);',3:4,1:2);
    c_eval('B?! = irf_tlim(B?,tints{!});',3:4,1:2);
    c_eval('B?!(:,2) = detrend(B?!(:,2));',3:4,1:2);
    %h=irf_plot({B3,B32trend,B32,B32filt})
    %irf_zoom(h,'x',tints{2})
    %irf_zoom(h,'y')
else
    c_eval('B?! = irf_add(1,staBpar?,-1,irf_filt(staBpar?,fmax,fmin,5));',3:4,1:2);
    %c_eval('B?! = irf_tlim(B?!,tints{!});',3:4,1:2);
    %c_eval('B?! = irf_add(1,B?!,-1,irf_filt(B?!,fmax,fmin));',3:4,1:2);
end

%
% simulation
toLoad = '20140802T231112-2007-08-31-500V_1600eV_9_5-5074743'; 
loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/'; 
doB = 0;

if doB    
    doOneR = 1;
    for sim = [0]; % fixed simulation
        ExB.run_V2B;
        switch sim
            case 0; modBz_rz = squeeze(mean(modBz,2))*1e9; BZ = modBz_rz; % nT
            case 1; simBz_rz = squeeze(mean(simBz,2))*1e9; BZ = simBz_rz;% nT
            case 2; fixBz_rz = squeeze(mean(fixBz,2))*1e9; BZ = fixBz_rz;% nT    
        end
    end
else 
    %load([loadPath toLoad],'modBz','simBz','zg','rg');
    load fixB_17km
    BZ = fixBz_rz;
end

load /Users/Cecilia/Data/BM/20070831/tcenter_2only.mat
c_eval('BB?! = [tocolumn(linspace(tints{!}(1),tints{!}(2),10000)) zeros(10000,1)];',3:4,1:2);

%c_eval('BB?! = [tocolumn(linspace(tints{!}(1),tints{!}(2),10000)) zeros(10000,1)];',3:4,1:2);
%BZ = modBz_rz(9,:);
%BZ = fixBz_rz;
distance = zg(2:end) - 0.5*diff(zg(1:2));

c_eval('timeseries?! = ts?(!)+distance/v!;',3:4,1:2);
c_eval('BB?! = [tocolumn(timeseries?!) tocolumn(BZ)];',3:4,1:2);
c_eval('BB?! = cn.zeropad(BB?!,tints{!});',3:4,1:2);
c_eval('BB?! = irf_filt(BB?!,fmin,fmax);',3:4,1:2);

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
nPanels = 8;
h = irf_plot(nPanels);
ww = 0.3; hh = 0.18; x1 = 0.19; x2 = 0.50;
yb = 0.1;
set(h(1),'position',[x1 yb+3*hh ww hh])
set(h(2),'position',[x1 yb+2*hh ww hh])
set(h(3),'position',[x1 yb+hh   ww hh])
set(h(4),'position',[x1 yb      ww hh])
set(h(5),'position',[x2 yb+3*hh ww hh])
set(h(6),'position',[x2 yb+2*hh ww hh])
set(h(7),'position',[x2 yb+hh   ww hh])
set(h(8),'position',[x2 yb      ww hh])
%
isub = 1;
%
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
    irf_plot(hca,Phi31,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Phi41,'color','b')
    irf_plot(hca,Phishift41,'b--')
    ylabel(hca,'\phi_{||} [V]')
    irf_zoom(hca,'y')
end
if 1 % magnetic field    
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,{B31,B41},'comp');
    ylabel(hca,'\delta B_{||} [nT]')   
end

% Second time intervalx
if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,E32,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,E42,'color','b')
    irf_plot(hca,Eshift42,'b--')
    ylabel(hca,'E_{||} [mV/m]')
    irf_zoom(hca,'y')
end
if 1 % Perpendicular electric field
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,Eper3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Eper4,'color','b')
    irf_plot(hca,Epershift42,'b--')
    ylabel(hca,'E_{\perp} [mV/m]')
    irf_zoom(hca,'y')
end
if 1 % Parallel electrostatic potential
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,Phi32,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Phi42,'color','b')
    irf_plot(hca,Phishift42,'b--')
    ylabel(hca,'\phi [V]')
    irf_zoom(hca,'y')
    
end
if 0 % magnetic field
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,BparAC3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,BparAC4,'color','b'); hold off;
    ylabel(hca,'\delta B_{||} [nT]')
end
if 1 % magnetic field    
    hca = h(isub); isub = isub + 1;  
    set(hca,'colororder',[0.2 0.5 0.2;0 0 0.9;1 0.8 0; 1 0 0.5])
    irf_plot(hca,{B32,B42,BB32,BB42},'comp');    
    hold(hca,'on')
    set(hca,'colororder',[0.2 0.5 0.2;0 0 0.9])
    irf_plot(hca,{BB32,BB42},'comp','linestyle',':'); 
    ylabel(hca,'\delta B_{||} [nT]');
    hold(hca,'off')
end
%

irf_zoom(h(1:4),'x',tints{1});
irf_zoom(h(5:8),'x',tints{2});
set(h([1 5]),'ylim',[-25 25])
set(h([2 6]),'ylim',[-4 19])
set(h([3 7]),'ylim',[-40 390])%set(h([3 6]),'ylim',[-70 240])
%irf_plot_axis_align;
irf_timeaxis(hca,'usefig');
for ff = 5:8
    set(h(ff),'yticklabel',[])
    ylabel(h(ff),' ')
end
tlabels1 = get(h(4),'xticklabel'); tlabels1{1} = ' '; tlabels1{11} = ' '; tlabels1{end} = ' ';
tlabels2 = get(h(8),'xticklabel'); tlabels2{1} = ' ';
set(h(8),'xticklabel',tlabels2)
set(h(4),'xticklabel',tlabels1)

if 1 % Ytick and Yticklabel ePhi/Te
    %%
    %set(gca,'Box','off');
    hca = h(7);
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
    ylabel(ax2,'e\phi/k_B T_{e}')
end
% common axis for b-field
Blim1 = get(h(4),'ylim');
Blim2 = get(h(8),'ylim');
blims = [min([Blim1 Blim2]) max([Blim1 Blim2])];
set(h(4),'ylim',blims);
set(h(8),'ylim',blims);

% Labeling
abcde={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.04,0.92],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:4%nPanels
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})    
    %set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    %irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end
for k=5:7%nPanels
    %irf_legend(h(k),x[abcde(k)],legloc{1},'color',colors{1})
    
    set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end
set(h(8),'colororder',[0.2 0.5 0.2;0 0 0.9;1 0.8 0; 1 0 0.5])
irf_legend(h(8),{'C3','C4','C3_{an}','C4_{an}'},[1 0.95])
%irf_legend(h(8),{'C3','C4','C3_{sim}','C4_{sim}'},[0.98 0.9])
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
    ax1=h(5);
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

if 0
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

end