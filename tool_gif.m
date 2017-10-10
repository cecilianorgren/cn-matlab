function [A,im,map] = tool_gif(x,y,z,correlation,phiE,ufEk,ufEn,Bz,title_str)

fig=figure(83);
set(gcf,'color','white'); % white background for figures (default is grey)
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
n_frames=max(size(correlation));
n_ims=100;
set(fig,'color','white'); 
set(fig,'position',[560   531   886   395])
h(1)=axes('position',[0.070    0.640    0.6750    0.270]);
h(3)=axes('position',[0.070    0.370    0.6750    0.270]);
h(4)=axes('position',[0.070    0.100    0.6750    0.270]);
h(2)=axes('position',[0.800    0.370    0.1250    0.2150]);

tint=[Bz(1,1) Bz(end,1)];
index0=find(correlation(:,1)==max(correlation(:,1)));
ind=0;

for k=fix(linspace(1,n_frames,n_ims));
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
    set(ax2,'YAxisLocation','right');
    set(ax2,'ytick',get(ax1,'ytick'))
    set(ax2,'yticklabel',num2str(get(ax1,'ytick')'*max(abs(Bz(:,2))),'%.2f'))
    set(ax2,'Color','none'); % color of axis
    set(ax2,'XColor','r','YColor','r'); % color of axis lines and numbers
    set(ax2,'xticklabel',[]);
    set(ax1,'xticklabel',[]);
    end
    title(h(1),title_str)
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
    %quiver3(h(2),0,0,0,n_hat(1),n_hat(2),n_hat(3),'g')
    plot3(h(2),x(:,1),x(:,2),x(:,3),'b')
    axis(h(2),'equal')
    set(h(2),'xlim',1.1*[-1 1],'ylim',1.1*[-1 1],'zlim',1.1*[-1 1])
    hold(h(2),'off')
    view(h(2),irf_cross(x(1,:),y(1,:)))
    corr_str=['corr = ',num2str(correlation(k),'%.003f')];
    b_str=['B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f')];
    vec_str_x=['k=[',num2str(x(k,1),'%0.1f'),' ',...
                   num2str(x(k,2),'%0.1f'),' ',...
                   num2str(x(k,3),'%0.1f'),']'];               
    vec_str_y=['n=[',num2str(y(k,1),'%0.1f'),' ',...
                   num2str(y(k,2),'%0.1f'),' ',...
                   num2str(y(k,3),'%0.1f'),']'];
    vec_str_z=['B=[',num2str(z(k,1),'%0.1f'),' ',...
                   num2str(z(k,2),'%0.1f'),' ',...
                   num2str(z(k,3),'%0.1f'),']'];
    title_right_str={b_str,vec_str_x,vec_str_y,vec_str_z,corr_str};
    %title(h(2),['B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f'),'\newline', corr_str]) %'\phi_{max}=',num2str(max(max(abs(phiE(:,2:end)))),'%.2f'),
    title(h(2),title_right_str)
    grid(h(2),'off')
    
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