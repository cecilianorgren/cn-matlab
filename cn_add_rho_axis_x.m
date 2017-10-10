set(gca,'box','off')
ax1=gca;
axesPosition=get(ax1,'Position');
xticks=get(gca,'xtick'); % in km
yticks=get(gca,'ytick'); % in km


yticks=yticks-yticks(1);
xticks=xticks-xticks(1);
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
ydiff=ylim(2)-ylim(1);
xdiff=xlim(2)-xlim(1);
n_yticks=size(yticks,2);
ytickstep=ceil(log10(ydiff/r_e));
xtickstep=ceil(log10(xdiff/r_e));

steps=fix(ydiff/r_e/ytickstep)
xsteps=fix(xdiff/r_e/xtickstep)
%
for k=1:xsteps
    xticks(k)=xlim(1)+2*k*r_e;
    xstr=num2str(2*k,'%.0f');
    xticklabels(k)={xstr};
end

axre=axes('Position',axesPosition,'ylimmode','auto','xlimmode','manual',...
        'yaxislocation','right',...
        'color','none','xaxislocation','top','xtick',xticks,...
        'xticklabel',xticklabels,'ytick',[],'yticklabel',[],...
        'Box','off','GridLineStyle','none',...
        'ylimmode','manual','ylim',ylim,'xlim',xlim);
xlabel(axre,'\rho_e')
linkaxes([ax1 axre])

if 0 % Ytick and Yticklabel B nT
    %set(gca,'Box','off');
    axesPosition=get(hcp,'Position');
    yticks0=get(hcp,'ytick');
    yticklabels_num=yticks0'/(6.2*1.6*0.1/Teav/scaling);
    yticklabels_str=num2str(yticklabels_num,'%.2f');
    ylim0=get(hcp,'ylim');
    xlim0=get(hcp,'xlim');
    ylim=ylim0/(6.2*1.6*0.1/Teav/scaling);
    axB=axes('Position',axesPosition,'ylimmode','manual','yaxislocation','right',...
        'color','none','xaxislocation','top','xtick',xticks0,...
        'xticklabel',[],'ytick',yticks0,'yticklabel',yticklabels_str,...
        'Box','off','GridLineStyle','none',...
        'ylimmode','manual','ylim',ylim0,'xlim',xlim0);
    ylabel(axB,'\delta B_z   [nT]')
end

