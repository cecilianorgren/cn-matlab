function [bB3 bB4] = ebplot(E3,E4,B3,B4,R3,R4,v,propdir,sE,sB)  
% +-y is propagation direction
% +-x is boundary normal direction

%% Get all E, B and R in same length
B3=irf_resamp(B3,E3);
B4=irf_resamp(B4,E3);
E4=irf_resamp(E4,E3);

% Center spacecraft positions
R3 = cn.mean(R3,1);
R4 = cn.mean(R4,1);
R0 = (R4 + R3)/2;
R3 = (R3 - R0); % km
R4 = (R4 - R0); % km

%% Spacecraft positions
dt = E3(2,1)-E3(1,1);
dy = v*dt; % km/s*s = km
% Now that we have dy (dt) we can take away frst column time series
B3 = B3(:,2:4);
B4 = B4(:,2:4);
E3 = E3(:,2:4);
E4 = E4(:,2:4);


fig = figure('name','Topology: E perp');
ax = plot(R3(1),R3(2),'go',R4(1),R4(2),'bv'); hold on;
xlabel('boundayr normal [km]'); ylabel('propagation direction [km]')


%% Plot Eperp and deltaB
red  =  [1 0 0];
blue =  [0 0.5 1];
green = [0 0.8 0];
line  = [0.7 0.7 0.7];

se = sE/max(E3(:,2));
sb = sB/max(B3(:,2));
%x0 = (Px3+Px4)/2;

%%
msize=5;

for l=1:1:size(E3,1) 
    if B3(l,3) > 0
        plot(R3(1),R3(2)+dy*(l-1),'ro','markersize',msize);
    else
        plot(R3(1),R3(2)+dy*(l-1),'rx','markersize',msize);
    end
    if B4(l,3) > 0
        plot(R4(1),R4(2)+dy*(l-1),'ro','markersize',msize);
    else
        plot(R4(1),R4(2)+dy*(l-1),'rx','markersize',msize);
    end
    %%
    if 1
    %plot(Pos3c(1),Pos3c(2),'go',Pos4c(2),Pos4c(1),'bv'); hold on;
    plot([R3(1) R4(1)],[R3(2) R4(2)]+dy*(l-1),'color',line,'linestyle',':')
    quiver(R4(1),R4(2)+dy*(l-1),B4(l,1)*sb,B4(l,2)*sb,0,'color',red)
    quiver(R3(1),R3(2)+dy*(l-1),B3(l,1)*sb,B3(l,2)*sb,0,'color',red)
    quiver(R4(1),R4(2)+dy*(l-1),E4(l,1)*se,E4(l,2)*se,0,'color',blue)
    quiver(R3(1),R3(2)+dy*(l-1),E3(l,1)*se,E3(l,2)*se,0,'color',green)
    
    %Px3=Px3+dx;%-dx;
    %Py3=Py3;
    %Px4=Px4+dx;
    %Py4=Py4;
    %pause(0.01)
    end
end
if 0
x1=(Px3+Px4-2*dx)/2;
if 0 % Add length scale for the E-arrows
    onekm=1;
    text(10,10,'10 km corresponds to ')
end

%% Adding scaling info
% Add an arrow that is 50 mV/m
% length 1 is scale*max(dE)
quiver(335,15,50*s,0,0,'color',[0 0 1]);
text(353,15,'50 mV/m')


%% Adding info to figure
time1=[num2str(t1(1)),'-',num2str(t1(2)),'-',num2str(t1(3)),...
    ' ',num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)),...
    '.',num2str(t1(7))];
time2=[num2str(t2(1)),'-',num2str(t2(2)),'-',num2str(t2(3)),...
    ' ',num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)),...
    '.',num2str(t2(7))];

if 1
    ax1=h(1);
    xlim=get(ax1,'Xlim');
    xticks0=get(ax1,'xtick');
    % find middle tick
    tick0=round(numel(xticks0)/2);
    % put each 3rd tick
    ticks = [flipdim(tick0:-2:0,2) 14:2:round(numel(xticks0))]
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
if 1
    ax1=gca;
    c4=text(435,-10,'C4');
    c3=text(460,10,'C3');
    xlabel(ax1,'x   [ km ]'); %xlabel(ax2,'x   [ \rho_e ]'); 
    ylabel(ax1,'y   [ km ]');
    %linkaxes([ax1 ax2],'x')
end

if 1
    %axis equal
    set(gca,'ylim',[-20 20])
    set(gca,'xlim',[330 490])
end
if 1
    xtmid=(x0+x1)/2;
    xt=get(ax1,'xtick');
    xt=xt+find((xt-xtmid)<0,1,'last');
    xtl={num2str((xt-xtmid)','%3.0f')};
    set(ax1,'xtick',xt,'xticklabel',xtl)
end
    
%title(['Wave perpendicular electric and parallel magnetic field']);

%legend('C3','C4');
set(gcf,'PaperPositionMode','auto'); 
ax1=gca;
if 1
cn_add_rho_axis_x
axis(ax1,'equal')
axis 'equal'
set(gca,'ylim',[-20 20])
end
%set(gca,'xlim',[330 490])
tit=title(gca,['\delta E_{\perp} and \delta B_{||}']);
set(gcf,'PaperPositionMode','auto');
%titpos=get(tit,'position');
%titpos(2)=titpos(2)*1.0
%set(tit,'position',titpos)
end
