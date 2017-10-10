% tool_boulder
% script to illustrate what the integrated potential is, depending on 
% assumed propagation direction, with comparison to wave magnetic field

if ~exist('gsmE3','var'); load matlabE; end
if ~exist('gsmB3','var'); load matlabB; end

% normal from MVA
n_hat=-[0.0856    0.9293   -0.3594];
n_hat=-[0.0856    -0.52   0.85];
sc=4;

for interval=10;
for f_filt=[1];
% time interval
switch interval
    case 1; t1= [2007 08 31 10 19 05.50]; t2=[2007 08 31 10 19 05.90]; % from article
    case 2; t1= [2007 08 31 10 19 07.20]; t2=[2007 08 31 10 19 07.50];
    case 3; t1= [2007 08 31 10 19 06.15]; t2=[2007 08 31 10 19 06.40];
    case 4; t1= [2007 08 31 10 19 06.90]; t2=[2007 08 31 10 19 07.50];
    case 5; t1= [2007 08 31 10 19 10.12]; t2=[2007 08 31 10 19 11.10]; 
    case 6; t1= [2007 08 31 10 19 04.50]; t2=[2007 08 31 10 19 05.10];
    case 7; t1= [2007 08 31 10 19 03.00]; t2=[2007 08 31 10 19 03.30]; % other?
    case 8; t1= [2007 08 31 10 19 07.50]; t2=[2007 08 31 10 19 08.00]; % other?
    case 9; t1= [2007 08 31 10 19 10.50]; t2=[2007 08 31 10 19 11.50]; % other?
    case 10; t1=[2007 08 31 10 18 42.45]; t2=[2007 08 31 10 18 43.55]; % dl
    %case 11; t1=[2007 08 31 10 39 28.00]; t2=[2007 08 31 10 39 30.00]; % other? doesnt go with tool_direction?
    case 11; t1=[2007 08 31 10 19 10.60]; t2=[2007 08 31 10 19 11.10]; % exist above
    case 12; t1=[2007 08 31 10 19 07.20]; t2=[2007 08 31 10 19 07.50]; % other?
    case 13; t1=[2007 08 31 10 18 45.10]; t2=[2007 08 31 10 18 45.20]; % other? perp?
end
tint=[toepoch(t1) toepoch(t2)];

% direction
k_direction=[1 0 0];

% get correlation, phiE, Bz and directions
n_tries=300;

c_eval('[x y z corr phiE Bz Ek En ufEn ufEk]=tool_direction(cn_toepoch(t1,t2,gsmB?),cn_toepoch(t1,t2,gsmE?),f_filt,n_tries);',sc);
index0=find(corr(:,1)==max(corr(:,1)));

%% make figures
clear im A f
n_ims=100;
fig=figure(83);
set(fig,'color','white'); 
setupfigure
set(fig,'position',[560   531   886   395])
h(1)=axes('position',[0.070    0.7100    0.6750    0.270]);
h(3)=axes('position',[0.070    0.4100    0.6750    0.270]);
h(4)=axes('position',[0.070    0.1100    0.6750    0.270]);
h(2)=axes('position',[0.800    0.1100    0.1250    0.2150]);
%%
ind=0;
for k=fix(linspace(1,n_tries,n_ims));
    ind=ind+1;   
    % potential match plot
    irf_plot(h(1),{[phiE(:,1) phiE(:,k+1)./repmat(max(max(abs(phiE(:,k+1)))),size(phiE,1),1)],...
        [Bz(:,1) Bz(:,2)./repmat(max(abs(Bz(:,2))),size(Bz,1),1)]},'comp'); 
    set(h(1),'ylim',[-1.1 1.1]);
    irf_legend(h(1),{'\phi_E','\phi_B'},[0.02 0.9]);
    if 0 % double axes with phi and B
    %ax1=h(1);
    %ax2 = axes('Position',get(ax1,'Position'));
    set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required    
    set(ax2,    'YAxisLocation','right');
    set(ax2,'ytick',get(ax1,'ytick'))
    set(ax2,'yticklabel',num2str(get(ax1,'ytick')'*max(abs(Bz(:,2))),'%.2f'))
    set(ax2,'Color','none'); % color of axis
    set(ax2,'XColor','r','YColor','r'); % color of axis lines and numbers
    set(ax2,'xticklabel',[]);
    set(ax1,'xticklabel',[]);
    end
    %title(h(1),['filtered at ',num2str(f_filt,'%.0f'),' Hz'])
    %irf_zoom(h(1),'x',tint); 
    grid(h(1),'off'); hold(h(1),'off')
    % electric field
    irf_plot(h(3),ufEk(:,[1 k+1])); ylabel(h(3),'E_k'); hold(h(3),'off')
    
    ylimk=[min(min(ufEk(:,2:end))) max(max(ufEk(:,2:end)))]; set(h(3),'ylim',ylimk);
    irf_plot(h(4),ufEn(:,[1 k+1])); ylabel(h(4),'E_n'); hold(h(4),'off')
    ylimn=[min(min(ufEn(:,2:end))) max(max(ufEn(:,2:end)))]; set(h(4),'ylim',ylimk);
    irf_zoom(h([1 3:4]),'x',tint);
    grid(h(3),'off');grid(h(4),'off')
    % direction plot    
    quiver3(h(2),0,0,0,x(k,1),x(k,2),x(k,3)); 
    hold(h(2),'on')
    quiver3(h(2),0,0,0,x(index0,1),x(index0,2),x(index0,3),'r')
    quiver3(h(2),0,0,0,n_hat(1),n_hat(2),n_hat(3),'g')
    plot3(h(2),x(:,1),x(:,2),x(:,3),'b')
    axis(h(2),'equal')
    set(h(2),'xlim',1.1*[-1 1],'ylim',1.1*[-1 1],'zlim',1.1*[-1 1])
    hold(h(2),'off')
    view(h(2),z(1,:))
    corr_str=['corr = ',num2str(corr(k),'%.003f')];
    title(h(2),['C3 \newline filtered at ',num2str(f_filt,'%.0f'),' Hz \newline' 'B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f'),'\newline', corr_str]) %'\phi_{max}=',num2str(max(max(abs(phiE(:,2:end)))),'%.2f'),
    grid(h(2),'off')
    
    % collect frame
    
    f=getframe(fig);
    A(:,ind)=f;
    if k==1, % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,n_ims) = 0;
    else
        im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
    end
    %pause(0.2)
end

%movie(A,n_ims);
%eval(['print -dpng /Users/Cecilia/EH/Pics/ANI_',varname,xtra_vn,'_1.png'])
%eval('imwrite(im,map,''/Users/Cecilia/Konferenser&Skolor/BoulderPostCluster12/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_',num2str(f_filt),'Hz,''DelayTime'',0,''LoopCount'',inf)');
imwrite(im,map,['/Users/Cecilia/Konferenser&Skolor/BoulderPostCluster12/',datestr(fromepoch(tint(1)),'HHMMSSFFF'),'_',num2str(f_filt),'Hz_PBE_',num2str(sc),'.gif'],'DelayTime',0.01,'LoopCount',inf);
close gcf;
end
end