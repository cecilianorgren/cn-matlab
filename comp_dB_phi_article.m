mu0=4*pi*1e-7;              % Vs/Am
B0=B(3);                    % nT
e=1.6e-19;                  % C
Ne=peaNeav*1e6;                % m^-3
%newNeav=0.20;
%Ne=newNeav*1e6;
%%
scaling1=(Ne*mu0*e/B0)*1e18; % gives RHS in nT or LHS in volt
scaling2=(peaNeav/(B0*5));
scaling=scaling1;
c_eval('LHS? = [bB?(:,1) bB?(:,4)/Teav/scaling];',3:4)

% Taking away average of phi
c_eval('RHS?=[phi?(:,1) (phi?(:,2)-mean(phi?(:,2)))/Teav];',3:4)
%%
eval(['LHS=LHS',Cstr,';']);
eval(['RHS=RHS',Cstr,';']);
%%
    fig=figure('name','Plottilotti');

if 1 % initialize figure
    set(fig,'Position',[10 400 500 750])
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'defaultlinelinewidth',1.0);
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
    %p(2)=axes('position',[0.55 0.77 0.3 0.2]); % [x y dx dy]
    
    tah=0.13;
    tah4=0.12;
    h(1)=axes('position',[0.15 tah4+3*(tah+0.005) 0.7 tah]); % [x y dx dy]
    h(2)=axes('position',[0.15 tah4+2*(tah+0.005) 0.7 tah]); % [x y dx dy]
    h(3)=axes('position',[0.15 tah4+1*(tah+0.005) 0.7 tah]); % [x y dx dy]
    h(4)=axes('position',[0.15 tah4 0.7 tah]); % [x y dx dy]
    
    p(1)=axes('position',[0.15 tah4+4*(tah+0.005)+0.03 0.3*1.1 0.3*1.1]); % [x y dx dy]
    %axes(p(1),square)
    
    ud=get(fig,'userdata');
    ud.subplot_handles=h;
    set(fig,'userdata',ud);
    set(fig,'defaultLineLineWidth',1);
    set(fig,'color','white')
end
if 1 % add sc configuration
    c_eval('arB?=cn_toepoch(t1,t2,gseB?);',3:4);
    [b_av b_hat b_mag]=cn_hat(arB3,arB4);
    
    z=b_hat;
    x=cn_cross(cn_cross(z,n_hat),z); % close to BL-normal
    y=cn_cross(z,x); % close to drift direction within boundary layer
    arM=[x;y;z];
    
    c_eval('arVi?=cn_toepoch(t1,t2,codifVi?);',3:4);
    [vi_av vi_hat vi_mag]=cn_hat(arVi3,arVi4);
    
    c_eval('arExB?=cn_toepoch(t1,t2,gseExB?);',3:4);
    [ExB_av ExB_hat ExB_mag]=cn_hat(arExB3,arExB4);
    
    t12=fix((t1+t2)/2);
    c_eval('argsePos?=cn_toepoch(t12,gsePos?);',3:4);

    % Transforming to field aligned system        
    Pos3=M*argsePos3(2:4)';
    Pos4=M*argsePos4(2:4)';
    Pos0=(Pos3+Pos4)/2;
    c_eval('Pos?c=Pos?-Pos0;',3:4);

    boxside=max(abs([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]));
    plot(p(1),Pos3c(1),Pos3c(2),'go','markersize',13,'linewidth',1.5); hold on;
    plot(p(1),Pos4c(1),Pos4c(2),'bv','markersize',13,'linewidth',1.5); hold on;
    c3=text(Pos3c(1),Pos3c(2),'C3   ','color','k','HorizontalAlignment','right');
    c4=text(Pos4c(1),Pos4c(2),'    C4','color','k','HorizontalAlignment','left');

    normdist=abs(Pos4c(1)-Pos3c(1));
    zdist=Pos4c(3)-Pos3c(3);
    kdist=abs(Pos4c(2)-Pos3c(2));
    cn_mag(Pos3c-Pos4c);

    xlabel('x   [km]'); ylabel('y   [km]');       
    axis([-boxside boxside -boxside boxside]*1.2)
    axis square       
    % Add B_0
    bcol=[0.6 0 0];
    bdot  = plot(p(1),-2,-2,'color',bcol,'linestyle','o','markersize',2,'linewidth',1);
    bring = plot(p(1),-2,-2,'color',bcol,'linestyle','o','markersize',20,'linewidth',1.1);
    b_txt = text(-1.5,-2,'  B','color',bcol);
    % Boundary normal
    ncol=[1 0 0.7];
    q1=quiver(0,0,1,0,4,'color',ncol); hold on;  
    nstr(1)={'boundary'};
    nstr(2)={'normal'};
    bn_txt=text(2.5,-1.2,nstr,'color',ncol);
    % Vi_av
    mvi=(arM*vi_av')';
    mvi_hat=(arM*vi_hat')';
    V_par=mvi(3);
    V_per=sqrt(mvi(1)^2+mvi(2)^2);
    q2=quiver(0,0,mvi_hat(1),mvi_hat(2),2*V_per/Vith,'color',[0.5 0.8 0.1]); hold on; 
    vistr(1)={'ion velocity '};
    vistr(2)={'from CIS'};
    vi_txt=text(-3,3.7,vistr,'color',[0.5 0.8 0.1]);
    % ExB_av
    mExB_hat=(arM*ExB_hat')';
    mExB_av=(arM*ExB_av')';
    mExB_per=sqrt(mExB_av(1)^2+mExB_av(2)^2);
    mExB_par=mExB_av(3);
    q3=quiver(0,0,mExB_hat(1),mExB_hat(2),2*mExB_per/Vith,'color',[1 0.6 0],'linewidth',1.1); hold on; 
    exb_txt=text(2,3,'ExB-velocity','color',[1 0.6 0]);
    % Propagation direction
    bn_hat=arM*n_hat';
    kcol=[0.5 0 1];
    q4=quiver(0,0,0,975,2/Vith,'color',[0.5 0 1]); hold on;
    kstr(1)={'wave'};
    kstr(2)={'velocity'};
    v_txt=text(-3.5,1.5,kstr,'color',kcol);
       
    set(bn_txt,'FontSize',10)
    set(vi_txt,'FontSize',10)
    set(exb_txt,'FontSize',10)
    %set(b_txt,'FontSize',10)
    set(v_txt,'FontSize',10)
    set(bn_txt,'FontSize',10)
    
    set(q1,'LineWidth',1.1);
    set(q2,'LineWidth',1.1);
    set(q3,'LineWidth',1.1);
    set(q4,'LineWidth',1.1);
end
if 1 % add info to figure
    clear infostr;
    nf=0;
    nf=nf+1;infostr(nf,1)={['\rho_e = ',num2str(round(r_e),'%.0f'),' km           \rho_i = ',num2str(10*round(r_i/10),'%.0f'),' km']};    
    nf=nf+1;infostr(nf,1)={['f_{LH} = ',num2str(round(flh),'%.0f'),' Hz']};
    nf=nf+1;infostr(nf,1)={['v_{||,i} = ',num2str(10*round(V_par/10),'%.0f'),' km/s     v_{\perp,i} = ',num2str(10*round(V_per/10),'%.0f'),' km/s']};
    nf=nf+1;infostr(nf,1)={['v_{ExB} = ',num2str(10*round(mExB_per/10),'%.0f'),' km/s']};
    nf=nf+1;infostr(nf,1)={['v_{th,i} = ',num2str(10*round(Vith/10),'%.0f'),' km/s']};   
    nf=nf+1;infostr(nf)={['T_{e} = ',num2str(100*round(Teav/100),'%.0f'),' eV']};    
    nf=nf+1;infostr(nf)={['n_{e} = ',num2str(peaNeav,'%.2f'),' cm^{-3}']};
    nf=nf+1;infostr(nf)={['\Delta t = ',num2str(t0choice,'%.4f'),' s']};
    nf=nf+1;infostr(nf)={['v_{wave} = ',num2str(10*round(vchoice/10),'%.0f'),' km/s']};
    nf=nf+1;infostr(nf,1)={['\Delta z_{C3,C4} = ',num2str(zdist,'%.0f'),' km']};
    nf=nf+1;infostr(nf)={['x = [',num2str(arM(1,1),'%.2f'),'  ',num2str(arM(1,2),'%.2f'),'  ',num2str(arM(1,3),'%.2f'),']  GSE']};
    nf=nf+1;infostr(nf)={['y = [',num2str(arM(2,1),'%.2f'),'  ',num2str(arM(2,2),'%.2f'),'  ',num2str(arM(2,3),'%.2f'),']  GSE']}; 
    nf=nf+1;infostr(nf)={['z = [',num2str(arM(3,1),'%.2f'),'  ',num2str(arM(3,2),'%.2f'),'  ',num2str(arM(3,3),'%.2f'),']  GSE']};     
    
    ann=annotation('textbox',[0.55 tah4+5*(tah+0.005)+0.16 0.5 0.2],...
    'BackgroundColor',[1 1 1],'Color',[0 0 0],...
    'String',infostr,'FontSize',10,...           
    'FitBoxToText','off','LineStyle','none',...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom'); 
end
if 0 % add figure with info, use above instead 
    set(p(2),'defaultTextFontSize',10);
    infostr(1)={['\rho_e = ',num2str(r_e,'%.1f'),' km']};
    infostr(2)={['f_{lh} = ',num2str(flh,'%.1f'),' Hz']};
    infostr(3)={['\Delta z_{C3,C4} = ',num2str(zdist,'%.0f'),' km']};
    %infostr(4)={['v_{||,CIS} = ',num2str(V_par,'%.0f'),' km/s']};
    %infostr(5)={['v_{\perp,CIS} = ',num2str(V_per,'%.0f'),' km/s']};
    %infostr(6)={['v_{ExB} = ',num2str(mExB_per,'%.0f'),' km/s']};
    infostr(4)={['v_{th,i} = ',num2str(Vith,'%.0f'),' km/s']};
    infostr(5)={['v_{wave} = ',num2str(vchoice,'%.0f'),' km/s']};
    text(p(2),0,1,infostr,'color',[0 0 0]),...
    %'HorizontalAlignment','left',...
    %'VerticalAlignment','top'); 
    plot(p(2),1,1,'xr')
    set(p(2),'xcolor','black','ycolor','white')
end

np=4;
isub=1;

if 1 % Electric field
    lb=[0.5 0.8 1]; % light blue
    lb=[0.6 0.7 1]; % light blue
    %lb=[0.5 0.8 0.2]; % light blue
    gr=[0 0.8 0]; % green
    
    c_eval('E?shift=bE?; E?shift(:,1)=E?shift(:,1)+t0choice;',4);
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,bE3(:,[1 2]),'color',gr); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 2]),'b'); hold(hca,'on');
    irf_plot(hca,E4shift(:,[1 2]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 2]),'b'); hold(hca,'on');  
    irf_plot(hca,bE3(:,[1 2]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    exlab=ylabel(hca,'E_{x}[mV/m]');
    set(hca,'ColorOrder',[[0 0.8 0];[0 0 1];lb;[0 0 0]]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.1]); 
    set(hca,'ColorOrder',[0 0 0]);  
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,bE3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 3]),'b'); hold(hca,'on');
    irf_plot(hca,E4shift(:,[1 3]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 3]),'b'); hold(hca,'on');  
    irf_plot(hca,bE3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    eylab=ylabel(hca,'E_{y}[mV/m]');
    set(hca,'ColorOrder',[[0 0.8 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.1]); 
    set(hca,'ColorOrder',[0 0 0]); 
    %set(hca,'ylim',[-30 80])
end
if 1 % Potential and dB
    hca=h(isub); isub=isub+1;
    hcp=hca;

    irf_plot(hca,RHS,'color',[0.8 0.5 0.2]); hold(hca,'on');    
    irf_plot(hca,LHS,'color',[0.4 0.2 0.8]); hold(hca,'on');
    set(hca,'ColorOrder',[[0.8 0.5 0.2];[0.4 0.2 0.8];[0 0 0]]);
    irf_legend(hca,{'\phi_{\delta E}','\phi_{\delta B} and \delta B'},[0.02 0.1])
    philab=ylabel(hca,['\phi_{\delta E}/T_e ,  \phi_{\delta B}/T_e'])
    set(hca,'ColorOrder',[0 0 0]);
    irf_legend(hca,{['C',Cstr]},[0.92 0.9])
end
if 1 % Gradient length scale
    hca=h(isub); isub=isub+1;

    irf_plot(hca,Ln2,'color',[1 0 0]); hold(hca,'on');
    %irf_plot(hca,Ln3,'color',[0 1 0]); hold(hca,'on');
    % irf_plot(hca,RHS,'color',[0.4 0.2 0.8]); hold on;
    %set(hca,'ColorOrder',[[0.8 0.5 0.2];[0.4 0.2 0.8]]);
    %irf_legend(hca,{''},[0.02 0.05])
    lnlab=ylabel(hca,'L_n / \rho_i');
    irf_zoom(hca,'y',[0 10])
    set(hca,'ylim',[0 3])
    
        
    ymstr(1)={'yellow marking'};
    ymstr(2)={'corresponds to'};
    ymstr(3)={'k_{perp}\rho_e=1   '};
    
    irf_legend(hca,{ymstr},[0.9 0.95],'FontSize',10,'color','black')
end
%
%irf_zoom('x',[cn_toepoch(t1) cn_toepoch(t2)]);
irf_zoom(h(1:4),'x',[phi3(1,1) phi3(end,1)])

if 1 % Xtick and Xticklabel length
    xticks0=get(hca,'xtick');
    midindex=ceil(length(xticks0)/2);
    xticklabels=cell(1,length(xticks0));
    for l=1:length(xticks0); xticklabels{l}=' '; end 
    if mod(length(xticks0-1),4)==0; si=1; else si=2; end
    
    for l=si:2:length(xticks0)
        xticklabels{l}= num2str((xticks0(l)-xticks0(midindex))*vchoice,'%.0f');
    end
    set(h(1),'xaxislocation','top','color','white','xtick',xticks0,'xticklabel',xticklabels);    
    xlabel(h(1),'km');
end
if 1 % Ytick and Yticklabel B nT
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

xy=get(gcf,'Position');
% set(gcf,'Position',[xy(1) xy(2) xy(3)*1.1 xy(4)*1.2])

abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)'};
for k=1:np
    irf_legend(h(k),[abcde(k)],[0.02,0.92],'color',[0 0 0])
end
if 1 % Marking one k*rho=1
    lambda=2*pi*r_e;
    deltat=lambda/vchoice;
    ind=size(phi3,1)/2;
    tfirst=phi3(fix(ind),1);

    tsecond=tfirst-deltat;
    for k=1:np
        irf_pl_mark(h(k),[tfirst tsecond],[1 1 0.0])
    end
    if 0
        tinterval=gsmB3(end,1)-gsmB3(1,1);
        ind=1;
        for k=1:size(gsmB3(end,1),1)
                vec(ind,1)=[gsmB3(k,1) 10]; 
                ind=ind+1;
          
        end
        irf_plot(h(2),vec)
    end
end
for k=1:np % Take away grid
    grid(h(k),'off')
end
if 0 % set all labels on same x-place
    philabpos=get(philab,'position')
    
    exlabpos=get(exlab,'position');
    exlabpos(1)=philabpos(1);
    set(exlab,'position')
    
    eylabpos=get(eylab,'position');
    eylabpos(1)=philabpos(1);
    set(eylab,'position')
    
    lnlabpos=get(lnlab,'position');
    lnlabpos(1)=philabpos(1);
    set(lnlab,'position')
    
end


set(gcf,'PaperPositionMode','auto');
eval(['print -depsc2 ',t1str,'_',t2str,'_dB_phi_article',Cstr,'.eps']);