function out = irf_match_phibe_vis(varargin)
% IRF_MATCH_PHIBE_VIS Visualizes the match from other match functions.
%   Used together with irf_match_phibe_dir.m and/or irf_match_phibe_v.m to
%   illustrate the matching that is made.
%
%   output = IRF_MATCH_PHIBE_VIS(type,req_in1,req_in2,...,op_in1,...)
% 
%   Examples:
%       % Direction
%       gif_stuff_dir = IRF_MATCH_PHIBE_VIS('direction',x,y,z,corr_dir,intEdt,Bz,En,Ek); 
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);
%
%       % Velocity
%       i_n=50; % if more than one densitiy, choose one by specifying index
%       gif_stuff_v = IRF_MATCH_PHIBE_VIS('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
%       imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);
%
%       % Density vs velocity
%       figure; h=axes;
%       axis_handle = IRF_MATCH_PHIBE_VIS('velocity/density',h,n,v,corr_v);
%
%   See also IRF_MATCH_PHIBE_DIR, IRF_MATCH_PHIBE_V

%% Check if axis is given 
[ax,args,nargs] = axescheck(varargin{:});
type = args{1};
varargin = args(2:end);
   
nargin;

switch type 
case 'direction' % gif with normalized phi_B and intEdt for different propagation directions    
    % Read input
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    corr_dir=varargin{4};
    A1=varargin{5}; % intEdt
    A2=varargin{6}; % Bz
    B=varargin{7}; % could be Ek,dEk
    C=varargin{8}; % could be En,dEn    
    dEn=varargin{9}; % dEn    
    dEk=varargin{10}; % dEk
    if nargin>11       
        mva_l = varargin{11};
        mva_v = varargin{12};
        f_filt = varargin{13}; % filter frequency
        f_filt_str = num2str(f_filt,'%.1f');
        B0 = varargin{14}; % B0
    else f_filt_str = ' unknown';
    end
    if ~exist('str_title','var')
        str_title=[];
    end
    % Initialize figure
    fig=figure;
    set(gcf,'color','white'); % white background for figures (default is grey)
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
    set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
    
    n_frames=max(size(corr_dir));
    n_ims=100; % make hundred images
    
    set(fig,'color','white'); 
    set(fig,'position',[560   531   886   395])
    h(1)=axes('position',[0.070    0.640    0.6300    0.270]);
    h(3)=axes('position',[0.070    0.370    0.6300    0.270]);
    h(4)=axes('position',[0.070    0.100    0.6300    0.270]);
    h(2)=axes('position',[0.800    0.470    0.1250    0.2150]); % small direction plot
    h(5)=axes('position',[0.820    0.100    0.100    0.2150]); % small direction plot

    tint=[A2(1,1) A2(end,1)];
    index0=find(corr_dir(:,1)==max(corr_dir(:,1))); % mark highest correlation with red indicator    
    ind=0;
    
    % make DC electric field, to have v_ExB
    dcEk = irf_add(1,B,-1,dEk);
    dcEn = irf_add(1,C,-1,dEn);

    try % find indices closest to mvars
        % max variance
        vecdiff = x-repmat(mva_v(1,:),size(x,1),1);
        dist_max = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarmax = find(dist_max==min(dist_max));
        % intermediate variance
        vecdiff = x-repmat(mva_v(2,:),size(x,1),1);
        dist_inter = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarinter = find(dist_inter==min(dist_inter));
        % minimum variance
        vecdiff = x-repmat(mva_v(3,:),size(x,1),1);
        dist_min = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarmin = find(dist_min==min(dist_min));
        ind_mvar = [ind_mvarmax ind_mvarinter ind_mvarmin];
        dist_mvar = [dist_max(ind_mvarmax) dist_inter(ind_mvarinter) dist_min(ind_mvarmin)];
    catch
        ind_mvar = [];
    end
    for k=fix(linspace(1,n_frames,n_ims));
        ind=ind+1;  
        
        % normalized potential match plot
        irf_plot(h(1),{[A1(:,1) A1(:,k+1)./repmat(max(max(abs(A1(:,k+1)))),size(A1,1),1)],...
            [A2(:,1) A2(:,2)./repmat(max(abs(A2(:,2))),size(A2,1),1)]},'comp'); 
        set(h(1),'ylim',[-1.1 1.1]);
        ylabel(h(1),'\phi/\phi_{max}')
        irf_legend(h(1),{'\phi_E','\phi_B'},[0.02 0.9]);
     
        title(h(1),str_title)      
        grid(h(1),'off'); hold(h(1),'off')

        % electric field
        irf_plot(h(3),{B(:,[1 k+1]),dEk(:,[1 k+1]),dcEk(:,[1 k+1])},'comp'); ylabel(h(3),'E_k [mV/m]'); hold(h(3),'off')    
        ylimk=[min(min(B(:,2:end))) max(max(B(:,2:end)))]; set(h(3),'ylim',ylimk);
        irf_legend(h(3),{'E','\delta E','E-\delta E'},[0.02 0.9]);
        irf_plot(h(4),{C(:,[1 k+1]),dEn(:,[1 k+1]),dcEn(:,[1 k+1])},'comp'); ylabel(h(4),'E_n [mV/m]'); hold(h(4),'off')    
        %irf_plot(h(4),C(:,[1 k+1])); ylabel(h(4),'E_n'); hold(h(4),'off')
        ylimn=[min(min(C(:,2:end))) max(max(C(:,2:end)))]; set(h(4),'ylim',ylimk);
        irf_legend(h(4),{'E','\delta E','E-\delta E'},[0.02 0.9]);
        irf_zoom(h([1 3:4]),'x',tint);
        grid(h(3),'off');grid(h(4),'off')

        % direction plot    
        quiver3(h(2),0,0,0,x(k,1),x(k,2),x(k,3)); 
        hold(h(2),'on')
        quiver3(h(2),0,0,0,x(index0,1),x(index0,2),x(index0,3),'r') % higgest correlation
        %quiver3(h(2),0,0,0,n_hat(1),n_hat(2),n_hat(3),'g')
        plot3(h(2),x(:,1),x(:,2),x(:,3),'b')
        axis(h(2),'equal')
        set(h(2),'xlim',1.1*[-1 1],'ylim',1.1*[-1 1],'zlim',1.1*[-1 1])
        
        view(h(2),irf_cross(x(1,:),y(1,:)))
        corr_str=['corr = ',num2str(corr_dir(k),'%.003f')];
        f_str=['f_{filt} = ', f_filt_str ' Hz'];
        b_str=['\delta B_{max}=',num2str(max(abs(A2(:,2))),'%.2f'),' nT'];
        vec_str_x=['k=[',num2str(x(k,1),'%0.1f'),' ',...
                       num2str(x(k,2),'%0.1f'),' ',...
                       num2str(x(k,3),'%0.1f'),']'];               
        vec_str_y=['n=[',num2str(y(k,1),'%0.1f'),' ',...
                       num2str(y(k,2),'%0.1f'),' ',...
                       num2str(y(k,3),'%0.1f'),']'];
        vec_str_z=['B=[',num2str(z(k,1),'%0.1f'),' ',...
                       num2str(z(k,2),'%0.1f'),' ',...
                       num2str(z(k,3),'%0.1f'),']'];
        title_right_str={b_str,f_str,vec_str_x,vec_str_y,vec_str_z};
        %title(h(2),['B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f'),'\newline', corr_str]) %'\phi_{max}=',num2str(max(max(abs(phiE(:,2:end)))),'%.2f'),
        title(h(2),title_right_str)
        
        % add minvar directions
        try            
            hmv1=quiver3(h(2),0,0,0,mva_v(1,1),mva_v(2,1),mva_v(3,1));
            hmv2=quiver3(h(2),0,0,0,mva_v(1,2),mva_v(2,2),mva_v(3,2));
            hmv3=quiver3(h(2),0,0,0,mva_v(1,3),mva_v(2,3),mva_v(3,3));
            mvarleg = legend(h(2),[hmv1 hmv2 hmv3],'max','inter','min','location','northwest');
            mvarleg.Position(1)=mvarleg.Position(1)+0.11;         
        end
        
        hold(h(2),'off')
        grid(h(2),'off')
        
        try
            if k ==1 % only need to add axes once
            % Add ExB on right axis
            if 0 
                ax1=h(3); % Ek, ExBn
                ax2 = axes('Position',get(ax1,'Position'));
                set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
                set(ax2,    'YAxisLocation','right');
                set(ax2,'Color','none'); % color of axis
                for k=1:numel(ax1.YTick)
                    yticklabels{k,1}=num2str(ax1.YTick(k)*(-round(1/B0*1e3)*10)/10,'%.f');
                end   
                set(ax2,'YTick',ax1.YTick,'YTickLabel',yticklabels)
                
                ax2.YTickLabel = yticklabels;
                
                if 0
                    xlim=get(ax1,'Xlim');
                    xticks0=get(ax1,'xtick');
                    % find middle tick
                    tick0=round(numel(xticks0)/2);
                    % put each 2rd tick
                    tickstep=2;
                    ticks = [flipdim(tick0:-tickstep:1,2) (tick0+tickstep):tickstep:round(numel(xticks0))];
                    tticks = xticks0(ticks);
                    lticks = round(tticks*v-xticks0(1)*v);   
                    for k=1:numel(lticks)
                        xticklabels{k}=num2str(lticks(k));
                    end    
                    xticklabels{end} = 'km'; 
                    yticks=[];
                    yticklabels=[];
                    ax2=axes('Position',get(ax1,'Position'),'Box','off',...
                             'XAxisLocation','top','YAxisLocation','right',...
                             'Color','none','Ytick',yticks,...
                             'YTickLabel',yticklabels,...
                             'xlim',xlim,'xtick',tticks([1:2:end-1 end]),...
                             'XTickLabel',xticklabels([1:2:end-1 end])); 
                    set(h(1),'box','off'); % remove 'xtick' if xticks required
                end
            else
                ax1=h(3); % Ek, ExBn
                ax2 = axes('Position',get(ax1,'Position'));
                set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
                set(ax2,    'YAxisLocation','right');
                set(ax2,'Color','none'); % color of axis
                set(ax2,'YLim',ax1.YLim/B0*1e3)
                set(ax2,'YTick',ax1.YTick/B0*1e3)
                for kk=1:numel(ax1.YTick)
                    yticklabels{kk,1}=num2str(ax1.YTick(kk)*(-round(1/B0*1e3)*10)/10,'%.f');
                end  
                set(ax2,'YTickLabel',yticklabels)
                
                % add label
                ax2.YLabel.String = '(ExB)_n [km/s]';
                
                ax1=h(4); % En, ExBn
                ax2 = axes('Position',get(ax1,'Position'));
                set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
                set(ax2,    'YAxisLocation','right');
                set(ax2,'Color','none'); % color of axis
                set(ax2,'YLim',ax1.YLim/B0*1e3)
                set(ax2,'YTick',ax1.YTick/B0*1e3)
                for kk=1:numel(ax1.YTick)
                    yticklabels{kk,1}=num2str(ax1.YTick(kk)*(round(1/B0*1e3)*10)/10,'%.f');
                end  
                set(ax2,'YTickLabel',yticklabels)
                
                % add label
                ax2.YLabel.String = '(ExB)_k [km/s]';
            end
            end
            
        end
        % correlation plot
        plot(h(5),1:numel(corr_dir),corr_dir,k,corr_dir(k),'go',index0,corr_dir(index0),'rx')        
        if 0 
            hold(h(5),'on')
            plot(h(5),ind_mvar(1)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(1)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(1)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(2)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(2)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(3)*[1 1],[-1 1],'color',hmv3.Color)
            plot(h(5),ind_mvar(3)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv3.Color)
            plot(h(5),ind_mvar(3)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv3.Color)            
        end
        hold(h(5),'off')
            
        ylabel(h(5),'corr')
        title(h(5),corr_str)
        set(h(5),'ylim',[-1 1],'xlim',[1 numel(corr_dir)],'xticklabel',[])
        
        % collect frames    
        f=getframe(fig);
        A(:,ind)=f;
        if k==1, % initialize animated gif matrix
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,n_ims) = 0;
        else
            im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    vis.A=A;
    vis.im=im;
    vis.map=map;
case 'velocity' % gif with phi_B and phi_E for different v 
    % Required input: phi_E,phi_B,v        
    phi_E=varargin{1};
    phi_B=varargin{2};    
    v=varargin{3};
    if ~isempty(varargin{4}); str_n=num2str(varargin{4},'%.2f'); else str_n='?'; end
    nv=length(v);
    
    % Adjust title_str and add density and B0
    if ~exist('title_str','var') 
        title_str=[];
    end
    %title_str=[title_str,', ',num2str(B0,'%.f'),' nT, ',num2str(n,'%.2f'), 'cc'];
    ylims=[floor(min(phi_B(:,end))/100) ceil(max(phi_B(:,end))/100)]*110;
    ylims=1.5*[min(phi_B(:,end)) max(phi_B(:,end))];

    fig=figure('name','Velocity match');
    set(gcf,'color','white'); % white background for figures (default is grey)
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    ind=0;
    for k=1:nv        
        h=irf_plot({phi_E(:,[1 1+k]),phi_B},'comp');
        str_legend{2}=['\phi_B(n=',str_n,'cc)'];
        str_legend{1}=['\phi_E(v=',num2str(v(k),'%0.f'),'km/s)'];                
        irf_legend(h,str_legend,[0.02 0.95]);        
        irf_zoom(h,'x',[phi_E(1,1) phi_E(end,1)]);
        set(h,'ylim',ylims); grid off;
        ylabel(h,'Potential [V]')
        title(h,title_str)
        hold off;

         % collect frames    
        f=getframe(fig);
        A(:,k)=f;
        if k==1, % initialize animated gif matrix
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,k) = 0;
        else
            im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    pause(1)
    vis.A=A;
    vis.im=im;
    vis.map=map;  
case 'best' % gif with phi_B and phi_E for different v         
    % Required input: phi_E,phi_B,v,corr_dir       
    phi_E=varargin{1};
    phi_B=varargin{2};    
    v=varargin{3};   
    k=varargin{4};   
    co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
    if ~isempty(ax)
        h = ax;        
        set(h,'ColorOrder',co)
    else
        fig=figure('name','Velocity match');
        set(gcf,'color','white'); % white background for figures (default is grey)
        set(gcf,'defaultAxesFontSize',14);
        set(gcf,'defaultTextFontSize',14);  
        set(fig,'position',[560   531   700   220])    
        h = irf_plot(1);
        h.Position = [0.12 0.20 0.80 0.60];
    end
    irf_plot(h,{phi_E,phi_B},'comp')   
    set(h,'ColorOrder',co)
    grid(h,'off')
    h(1).YLabel.String = '\phi [V]';
    irf_legend(h,{'\phi_E','\phi_B'},[0.98, 0.95]);
    irf_legend(h,['v = [' num2str(k,'%.2f') ']\times' num2str(v,'%.f') ' km/s'],[0.02, 0.95]);    
    irf_zoom(h,'x',[phi_E(1,1) phi_E(end,1)])
    irf_zoom(h,'y')
    
    
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
        lticks = round(tticks*v-xticks0(1)*v);   
        for k=1:numel(lticks)
            xticklabels{k}=num2str(lticks(k));
        end    
        xticklabels{end} = 'km'; 
        yticks=[];
        yticklabels=[];
        ax2=axes('Position',get(ax1,'Position'),'Box','off',...
                 'XAxisLocation','top','YAxisLocation','right',...
                 'Color','none','Ytick',yticks,...
                 'YTickLabel',yticklabels,...
                 'xlim',xlim,'xtick',tticks([1:2:end-1 end]),...
                 'XTickLabel',xticklabels([1:2:end-1 end])); 
         set(h(1),'box','off'); % remove 'xtick' if xticks required
    end

    %title_str=[title_str,', ',num2str(B0,'%.f'),' nT, ',num2str(n,'%.2f'), 'cc'];
    %ylims=[floor(min(phi_B(:,end))/100) ceil(max(phi_B(:,end))/100)]*110;
    %ylims=1.5*[min(phi_B(:,end)) max(phi_B(:,end))];

    if nargin>5 
        corr_dir=varargin{5};
        prelpos = h.Position;
        h.Position = [0.12 prelpos(2) 0.55 prelpos(4)];        
        ax2.Position = [0.12 prelpos(2) 0.55 prelpos(4)];
        hc = axes;
        hc.Position = [0.85 [prelpos(2) 0.12 prelpos(4)]];
        index0=find(corr_dir(:,1)==max(corr_dir(:,1))); % mark highest correlation with red indicator    
        plot(hc,1:numel(corr_dir),corr_dir,index0,corr_dir(index0),'rx')        
        hc.XTickLabel = [];
        hc.XLim = [1 numel(corr_dir)];
        title(hc,['C_{max} = ' num2str(corr_dir(index0),'%.2f')])
        hc.YLabel.String = 'C';
    end
    if nargin>6 
        flim=varargin{6};
        flh=varargin{7};
        irf_legend(h,['f_{filt} = ' num2str(flim,'%.2f') '\times f_{LH} = ' num2str(flim*flh,'%.f') ' Hz'],[0.02, 0.09]);
        h.YLim = h.YLim*1.2;
    end 
    try; if nargin>8 % ratio between phi_B and Bz, for yaxis
        yratio=varargin{8};        
        %ax2.YTick = ax1.YTick/yratio;
        ax2.YLim = ax1.YLim/yratio;
        ax2.YLabel.String = 'B_{||} [nT]'; 
        ax2.YTickLabelMode = 'auto';
        ax2.YTickMode = 'auto';
    end; end 
    try; if nargin>9 
        re=varargin{9};        
        irf_legend(h,['\rho_e = ' num2str(re,'%.1f') 'km'],[0.98, 0.09]);
        h.YLim = h.YLim*1.2;
    end; end 
    vis = [];
case 'velocity/density' % 2D plot of correlation=f(v,n)
    [ax,args,nargs] = axescheck(varargin{:});
    if isempty(ax); ax=axes; end
    n=args{1};
    v=args{2};
    corr_v=args{3};
    pcolor(ax,v,n,log10(corr_v)); 
    shading(ax,'flat')
    title(ax,'Correlation')
    xlabel(ax,'Velocity')
    ylabel(ax,'Density')
    axc=colorbar('peer',ax);               
    ylabel(axc,'log10(sum((\phi_E-\phi_B)^2))')    
    vis.ax=ax;
    vis.axc=axc;
case 'magnetic field' % gif with normalized phi_B and intEdt for different propagation directions    
    % Read input
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    corr_dir=varargin{4};
    A1=varargin{5}; % intEdt
    A2=varargin{6}; % Bz
    B=varargin{7}; % could be Ek,dEk
    C=varargin{8}; % could be En,dEn    
    dEn=varargin{9}; % dEn    
    dEk=varargin{10}; % dEk
    if nargin>11       
        mva_l = varargin{11};
        mva_v = varargin{12};
        f_filt = varargin{13}; % filter frequency
        f_filt_str = num2str(f_filt,'%.1f');
    else f_filt_str = ' unknown';
    end
    if ~exist('str_title','var')
        str_title=[];
    end
    % Initialize figure
    fig=figure;
    set(gcf,'color','white'); % white background for figures (default is grey)
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
    set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
    
    n_frames=max(size(corr_dir));
    n_ims=100; % make hundred images
    
    set(fig,'color','white'); 
    set(fig,'position',[560   531   886   395])
    h(1)=axes('position',[0.070    0.640    0.6750    0.270]);
    h(3)=axes('position',[0.070    0.370    0.6750    0.270]);
    h(4)=axes('position',[0.070    0.100    0.6750    0.270]);
    h(2)=axes('position',[0.800    0.470    0.1250    0.2150]); % small direction plot
    h(5)=axes('position',[0.820    0.100    0.100    0.2150]); % small direction plot

    tint=[A2(1,1) A2(end,1)];
    index0=find(corr_dir(:,1)==max(corr_dir(:,1))); % mark highest correlation with red indicator    
    ind=0;
    
    % make DC electric field, to have v_ExB
    dcEk = irf_add(1,B,-1,dEk);
    dcEn = irf_add(1,C,-1,dEn);

    try % find indices closest to mvars
        % max variance
        vecdiff = x-repmat(mva_v(1,:),size(x,1),1);
        dist_max = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarmax = find(dist_max==min(dist_max));
        % intermediate variance
        vecdiff = x-repmat(mva_v(2,:),size(x,1),1);
        dist_inter = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarinter = find(dist_inter==min(dist_inter));
        % minimum variance
        vecdiff = x-repmat(mva_v(3,:),size(x,1),1);
        dist_min = sqrt(vecdiff(:,1).^2+vecdiff(:,2).^2+vecdiff(:,3).^2);
        ind_mvarmin = find(dist_min==min(dist_min));
        ind_mvar = [ind_mvarmax ind_mvarinter ind_mvarmin];
        dist_mvar = [dist_max(ind_mvarmax) dist_inter(ind_mvarinter) dist_min(ind_mvarmin)];
    catch
        ind_mvar = [];
    end
    for k=fix(linspace(1,n_frames,n_ims));
        ind=ind+1;  
        
        % normalized potential match plot
        irf_plot(h(1),{[A1(:,1) A1(:,k+1)./repmat(max(max(abs(A1(:,k+1)))),size(A1,1),1)],...
            [A2(:,1) A2(:,2)./repmat(max(abs(A2(:,2))),size(A2,1),1)]},'comp'); 
        set(h(1),'ylim',[-1.1 1.1]);
        ylabel(h(1),'\phi/\phi_{max}')
        irf_legend(h(1),{'\phi_E','\phi_B'},[0.02 0.9]);
     
        title(h(1),str_title)      
        grid(h(1),'off'); hold(h(1),'off')

        % electric field
        irf_plot(h(3),{B(:,[1 k+1]),dEk(:,[1 k+1]),dcEk(:,[1 k+1])},'comp'); ylabel(h(3),'E_k'); hold(h(3),'off')    
        ylimk=[min(min(B(:,2:end))) max(max(B(:,2:end)))]; set(h(3),'ylim',ylimk);
        irf_legend(h(3),{'E','\delta E','E-\delta E'},[0.02 0.9]);
        irf_plot(h(4),{C(:,[1 k+1]),dEn(:,[1 k+1]),dcEn(:,[1 k+1])},'comp'); ylabel(h(4),'E_n'); hold(h(4),'off')    
        %irf_plot(h(4),C(:,[1 k+1])); ylabel(h(4),'E_n'); hold(h(4),'off')
        ylimn=[min(min(C(:,2:end))) max(max(C(:,2:end)))]; set(h(4),'ylim',ylimk);
        irf_legend(h(4),{'E','\delta E','E-\delta E'},[0.02 0.9]);
        irf_zoom(h([1 3:4]),'x',tint);
        grid(h(3),'off');grid(h(4),'off')

        % direction plot    
        quiver3(h(2),0,0,0,x(k,1),x(k,2),x(k,3)); 
        hold(h(2),'on')
        quiver3(h(2),0,0,0,x(index0,1),x(index0,2),x(index0,3),'r') % higgest correlation
        %quiver3(h(2),0,0,0,n_hat(1),n_hat(2),n_hat(3),'g')
        plot3(h(2),x(:,1),x(:,2),x(:,3),'b')
        axis(h(2),'equal')
        set(h(2),'xlim',1.1*[-1 1],'ylim',1.1*[-1 1],'zlim',1.1*[-1 1])
        
        view(h(2),irf_cross(x(1,:),y(1,:)))
        corr_str=['corr = ',num2str(corr_dir(k),'%.003f')];
        f_str=['f_{filt} = ', f_filt_str ' Hz'];
        b_str=['\delta B_{max}=',num2str(max(abs(A2(:,2))),'%.2f'),' nT'];
        vec_str_x=['k=[',num2str(x(k,1),'%0.1f'),' ',...
                       num2str(x(k,2),'%0.1f'),' ',...
                       num2str(x(k,3),'%0.1f'),']'];               
        vec_str_y=['n=[',num2str(y(k,1),'%0.1f'),' ',...
                       num2str(y(k,2),'%0.1f'),' ',...
                       num2str(y(k,3),'%0.1f'),']'];
        vec_str_z=['B=[',num2str(z(k,1),'%0.1f'),' ',...
                       num2str(z(k,2),'%0.1f'),' ',...
                       num2str(z(k,3),'%0.1f'),']'];
        title_right_str={b_str,f_str,vec_str_x,vec_str_y,vec_str_z};
        %title(h(2),['B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f'),'\newline', corr_str]) %'\phi_{max}=',num2str(max(max(abs(phiE(:,2:end)))),'%.2f'),
        title(h(2),title_right_str)
        
        % add minvar directions
        try            
            hmv1=quiver3(h(2),0,0,0,mva_v(1,1),mva_v(2,1),mva_v(3,1));
            hmv2=quiver3(h(2),0,0,0,mva_v(1,2),mva_v(2,2),mva_v(3,2));
            hmv3=quiver3(h(2),0,0,0,mva_v(1,3),mva_v(2,3),mva_v(3,3));
            mvarleg = legend(h(2),[hmv1 hmv2 hmv3],'max','inter','min','location','northwest');
            mvarleg.Position(1)=mvarleg.Position(1)+0.11;         
        end
        
        hold(h(2),'off')
        grid(h(2),'off')

        % correlation plot
        plot(h(5),1:numel(corr_dir),corr_dir,k,corr_dir(k),'go',index0,corr_dir(index0),'rx')        
        if 0 
            hold(h(5),'on')
            plot(h(5),ind_mvar(1)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(1)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(1)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv1.Color)
            plot(h(5),ind_mvar(2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(2)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(2)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv2.Color)
            plot(h(5),ind_mvar(3)*[1 1],[-1 1],'color',hmv3.Color)
            plot(h(5),ind_mvar(3)-fix(n_frames/2)*[1 1],[-1 1],'color',hmv3.Color)
            plot(h(5),ind_mvar(3)+fix(n_frames/2)*[1 1],[-1 1],'color',hmv3.Color)            
        end
        hold(h(5),'off')
            
        ylabel(h(5),'corr')
        title(h(5),corr_str)
        set(h(5),'ylim',[-1 1],'xlim',[1 numel(corr_dir)],'xticklabel',[])
        
        % collect frames    
        f=getframe(fig);
        A(:,ind)=f;
        if k==1, % initialize animated gif matrix
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,n_ims) = 0;
        else
            im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    vis.A=A;
    vis.im=im;
    vis.map=map;

end

out=vis;