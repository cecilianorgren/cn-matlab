
h=irf_plot({E3,E4},'comp');

for kk=1:3; 
    grid(h(kk),'off'); 
    irf_legend(h(kk),{'C3','C4'},[0.02 0.95]);
    irf_legend(h(kk),{'GSE'},[0.02 0.75],'color',[0  0 0]);
end
irf_zoom(h,'y')
irf_zoom(h,'x',tint)

ylabel(h(1),'E_x [mV/m]')
ylabel(h(2),'E_y [mV/m]')
ylabel(h(3),'E_z [mV/m]')

if 1
    v1=450;
    v2=450;
    %
    % Left panels
    %
    ax1=h(1);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    tick0=2;
    % put each 2rd tick
    tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    tticks = xticks0(ticks);
    lticks = round(tticks*v1-1*xticks0(2)*v1);
    %lticks = round(tticks*v1-xticks0(tick0)*v1);
    
    %nticksl=xlim(1)/5/r_e;
    %nticksr=xlim(2)/5/r_e;
    %xticks=0:xlim(2)/nticksr:xlim(2);
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end    
    %xticklabels{1} = ' '; 
    %xticklabels{2} = ' ';
    xticklabels{end} = ' '; 
    xticklabels{end} = 'km'; 
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
    if 0
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
end