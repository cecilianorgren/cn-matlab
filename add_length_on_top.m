function ax2 = add_length_on_top(ax,v,tickstep)

if 0
    ax1=ax;
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    % put each 2rd tick
    %tickstep=2;
    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
    tticks = xticks0(ticks);
    lticks = round(tticks*v-xticks0(1)*v);   
    for k=1:numel(lticks)
        xticklabels{k}=num2str(lticks(k));
    end    
    %xticklabels{end} = 'km'; 
    yticks=[];
    yticklabels=[];
    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
             'XAxisLocation','top','YAxisLocation','right',...
             'Color','none','Ytick',yticks,...
             'YTickLabel',yticklabels,...
             'xlim',xlim,'xtick',tticks,...
             'XTickLabel',xticklabels);
    set(ax,'box','off'); % remove 'xtick' if xticks required

else
    ax1=ax;
    xticks0=get(ax1,'xtick');
    
    xlim=get(ax1,'Xlim');
    T = diff(xlim);
    L = T*v;
    dL = 10^ceil(log10(L))/10;
    Lticks = 0:(dL*tickstep):L;
    tticks = Lticks/v+xticks0(1);
          
    for k=1:numel(Lticks)
        xticklabels{k}=num2str(Lticks(k)-mean(Lticks));
    end    
    %xticklabels{end} = 'km'; 
    yticks=[];
    yticklabels=[];
    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
             'XAxisLocation','top','YAxisLocation','right',...
             'Color','none','Ytick',yticks,...
             'YTickLabel',yticklabels,...
             'xlim',xlim,'xtick',tticks,...
             'XTickLabel',xticklabels);
     set(ax,'box','off'); % remove 'xtick' if xticks required                       
end