%% eh_running
cd /Users/Cecilia/Data/BM/20070831
tint=eh_tint;
ntint=length(tint);
% prepare peace data and load E and B
load matlabE;
load matlabB;
c_eval('peace_e? = c_caa_distribution_data(''C?_CP_PEA_3DXPH_DEFlux'');',3:4)
% make E parallel
angle_lim=90;
c_eval('[E?,d?]=irf_edb(diE?,diB?,angle_lim,''Epar'');',3:4);
c_eval('E?proj = irf_dot(E?,irf_norm(diB?));',3:4)
c_eval('E?noproj = irf_dot([diE?(:,1:3) diE?(:,1)*0],irf_norm(diB?));',3:4)
%%
save_path = '/Users/Cecilia/Research/EH/2012-04-22/';
plot_type={'polar','cross-section'};

%for p=1:1; % both plot types
for k=1:length(tint);
    emin=100;
    emin_scale=emin;
    fig=figure(101);
    clf(fig)
    setupfigure
    set(fig,'Renderer','painters')
    set(fig,'color','white')
    set(fig,'PaperUnits','centimeters')
    xSize = 15; ySize = 15;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(fig,'PaperPosition',[xLeft yTop xSize ySize])
    pos=get(fig,'Position');
    set(fig,'Position',[pos(1) pos(2) xSize*60 ySize*60])
    ax1=axes('position',[0.15 0.6 0.3 0.32]); % [x y dx dy]
    ax2=axes('position',[0.15 0.15 0.3 0.32]); % [x y dx dy]
    ax3=axes('position',[0.55 0.52 0.35 0.45]); % [x y dx dy]
    ax4=axes('position',[0.55 0.3 0.3 0.1]); % [x y dx dy]
    ax5=axes('position',[0.55 0.15 0.3 0.1]); % [x y dx dy]
    setupfigure
%
    [~,cb]=c_caa_plot_distribution_function(ax1,'tint',tint{k},...
        plot_type{2},'emin_scale',emin_scale,'emin',emin,peace3);
    [~,cb]=c_caa_plot_distribution_function(ax2,'tint',tint{k},...
        plot_type{2},'emin_scale',emin_scale,'emin',emin,peace4);
    [~,cb]=c_caa_plot_distribution_function(ax3,'tint',tint{k},...
        plot_type{1},'emin_scale',emin_scale,'emin',emin,peace3,peace4);
    
    set(ax1,'xlim',get(ax2,'xlim'),'ylim',get(ax2,'ylim'))
    %set(ax1,'position',[0.15 0.62 0.7 0.32])
    %set(cb,'ylim',[-3 1])
    irf_plot(ax4,{irf_tlim(E3(:,[1 4]),tint{k}),irf_tlim(E4(:,[1 4]),tint{k})},'comp')
    irf_legend(ax4,{'C3','C4'},[0.02,0.90])
    grid(ax4,'off')
    irf_zoom(ax4,'x',tint{k})
    ylabel(ax4,'E_{||}')
    
    irf_plot(ax5,{[diE3(:,1),d3],[diE4(:,1), d4]},'comp')
    grid(ax5,'off')
    irf_legend(ax5,{'C3','C4'},[0.02,0.90])
    irf_zoom('x',tint{k}) % [diE3(1,1) diE3(end,1)]
    ylabel(ax5,'B-Spinplane')
    print('-dpng',[save_path,'pol_cs_e_ang_34_',num2str(k)])
end
%end

%% make similar plots but for all times, too find the strongest beam, 5s intervals
% find tmin and tmax from tint
tstart=toepoch([2007 08 31 10 15 00]);
tstop=toepoch([2007 08 31 10 20 00]);

tstep=4; % 4s = one spin
emin=100;
emin_scale=emin;
fn=figure(99);
setupfigure
t1=tstart;
t2=tstart+4;
count=1;
%%
while t1 < tstop
    for k=1:4;
        tt=[t1 t2]; t1=t1+tstep; t2=t2+tstep;
        ax(k)=subplot(2,2,k);
        [axx(k),cb(k)]=c_caa_plot_distribution_function(ax(k),'tint',tt,...
            'polar','emin_scale',emin_scale,'emin',emin,peace_e3,peace_e4);        
    end
    linkaxes(axx)
    set(axx,'clim',[5 7.5])
    print('-dpng',[save_path,'pol_34_allt_deflux_ax_',num2str(count)])
    clf(fn);
    count=count+1;
end