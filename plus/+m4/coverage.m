% assume we have functions x(t), y(t) in variables x, y, t, in [m] and [s]
RE = 6371*1e3; % m, Earth radius

r_aps = ([15 25 61]+1)*RE; % m, apogee
r_pers = ([3 3 13]+1)*RE; % m, perigee
t=0; x=[]; y=[];
dt = 120; % s, timestep, 2 min 
Ttot = 60*60*24*365; % one year

for year=[1 2 3]; 
    [tYear,xYear,yYear] = m4.orbit(r_pers(year),r_aps(year),dt,Ttot);
    t = [t tYear+t(end)];
    x = [x xYear];
    y = [y yYear];
end
t=t(2:end);
   %%     
        
        
        
% Set up grid    
rmax = ceil(max([x y])*1.1/RE)*RE; % m
binsize = 1*RE; % m
nbins = 2*rmax/binsize;
edges = linspace(-rmax,rmax,nbins+1);
centers = (edges(1:nbins)+binsize/2);
dt = diff(t(1:2));

if 1 % illustrate gridding
    [X,Y] = meshgrid(edges,edges);
    NT = histcn([x(:) y(:)], edges, edges);   
    TT = NT*dt;
    Ttot = sum(sum(TT));
    ndays = sum(sum(TT))/60/60/24;
    
    hca=subplot(1,3,1);
    hold(hca,'on');    
    mesh(hca,X/RE,Y/RE,X*0-1);
    plot(hca,x/RE,y/RE,'b')
    axis(hca,'equal')
    set(hca,'xlim',rmax/RE*[-1 1],'ylim',rmax/RE*[-1 1])
    title(hca,['Orbit during T_{tot} = ' num2str(ndays,'%.0f') ' days'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E')
    
    hca=subplot(1,3,2);
    hold(hca,'on');    
    mesh(hca,X/RE,Y/RE,X*0-1);
    plot(hca,x/RE,y/RE,'b',x/RE,y/RE,'bo')
    axis(hca,'equal')
    set(hca,'xlim',[10 12],'ylim',[10 12])
    title(hca,['time between two consecutive ''o'', dt = ' num2str(dt,'%.0f') ' s'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E')  
    
    hca=subplot(1,3,3);
    ax_pos=get(hca,'position');   
    [cX,cY] = meshgrid(centers,centers);
    %surf(hca,cX/RE,cY/RE,TT/60/60)
    surf(hca,X/RE,Y/RE,X*0,TT/60)
    view(hca,[0 0 1])
    axis(hca,'equal')
    set(hca,'xlim',rmax/RE*[-1 1],'ylim',rmax/RE*[-1 1])%,'clim',dt*[0 20])
    title(hca,['Orbit coverage, T_{tot} = ' num2str(ndays,'%.0f') ' days'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E') 
    shading(hca,'flat')
    
    cb = colorbar('peer',hca);
    ylabel(cb,['total time spent in bin [min]'])
    set(hca,'position',ax_pos);
    
    %%
    % change colormap so that 0 is black
    cmap=colormap('gray');
    
    if 1 % make background white        
        cmap=flipdim(cmap,1);
        cmap(2:3,:)=[];
    else
        cmap(2:10,:)=[];
    end
    %cmap(1,:)=[0 0 0];
    colormap(cmap);
    
    
    %cmap((end-1):(end-3),:)=[];
    
    %colormap(cmap);
end
