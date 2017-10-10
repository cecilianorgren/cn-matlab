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
figure('name','dB and phi');
np=3;
h=irf_plot(np);
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
    
    ylabel(hca,'E_{x}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,bE3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 3]),'b'); hold(hca,'on');
    irf_plot(hca,E4shift(:,[1 3]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,bE4(:,[1 3]),'b'); hold(hca,'on');  
    irf_plot(hca,bE3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{y}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]); 
end
if 1 % Potential and dB
    hca=h(isub); isub=isub+1;

    irf_plot(hca,LHS,'color',[0.8 0.5 0.2]); hold on;    
    irf_plot(hca,RHS,'color',[0.4 0.2 0.8]); hold on;
    irf_zoom('x',[cn_toepoch(t1) cn_toepoch(t2)]);
    set(hca,'ColorOrder',[[0.8 0.5 0.2];[0.4 0.2 0.8]]);
    irf_legend(hca,{'phi_{\delta B} and \delta B','phi_{\delta E}'},[0.02 0.05])
    ylabel(['     phi/T_e  ,    phi_{\delta B}/T_e    C',Cstr])
end

%title(['Electrostatic potential and magnetic field fluctuations  \newline (n_e=',num2str(peaNeav,'%.3f'),' cm^{-3} B_0=',num2str(B0,'%.1f'),' nT)'])
irf_zoom(h,'x',[phi3(1,1) phi3(end,1)])

if 1 % Xtick and Xticklabel length
    xticks0=get(hca,'xtick');
    midindex=ceil(length(xticks0)/2);
    xticklabels=cell(1,length(xticks0));
    for l=1:length(xticks0); xticklabels{l}=' '; end 
    if mod(length(xticks0-1),4)==0; si=1; else si=2; end
    
    for l=si:2:length(xticks0)
        xticklabels{l}= num2str((xticks0(l)-xticks0(midindex))*v(f,g),'%.0f');
    end
    set(h(1),'xaxislocation','top','color','white','xtick',xticks0,'xticklabel',xticklabels);    
    xlabel(h(1),'km');
end
if 1 % Ytick and Yticklabel B nT
    set(gca,'Box','off');
    axesPosition=get(h(3),'Position');
    yticks0=get(h(3),'ytick');
    yticklabels_num=yticks0'/(6.2*1.6*0.1/Teav/scaling);
    yticklabels_str=num2str(yticklabels_num,'%.2f');
    ylim0=get(hca,'ylim');
    xlim0=get(hca,'xlim');
    ylim=ylim0/(6.2*1.6*0.1/Teav/scaling);
    axB=axes('Position',axesPosition,'ylimmode','manual','yaxislocation','right',...
        'color','none','xaxislocation','top','xtick',xticks0,...
        'xticklabel',[],'ytick',yticks0,'yticklabel',yticklabels_str,...
        'Box','off',...
        'ylimmode','manual','ylim',ylim0,'xlim',xlim0);
    ylabel(axB,'\delta B_z   [nT]')
end

title(h(1),['\Delta t=',...
    num2str(t0choice,'%.4f'),' s   v=',...
    num2str(v(f,g),'%.0f'),' km/s   T_e=',...
    num2str(Teav,'%.0f'),' eV   B_0=',...
    num2str(B0,'%.0f'),' nT   n_e=',...
    num2str(peaNeav,'%.2f'),' cm^{-3}']);
%'Filtered and shited electric field, electrostatic potential and filtered magnetic field')
xy=get(gcf,'Position');
set(gcf,'Position',[xy(1) xy(2) xy(3)*1.1 xy(4)*1.2])

abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)'};
for k=1:np
    irf_legend(h(k),[abcde(k)],[0.02,0.92],'color',[0 0 0])
end
if 1 % Marking one k*rho=1
    lambda=2*pi*r_e;
    deltat=lambda/vchoice;
    tfirst=phi3(fix(end/2),1);

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
    %markvec(:,1)=[tfirst;tsecond];
    %markvec(:,2)=[60 60];
    %irf_plot(h(2),markvec,'+-')
end


set(gcf,'PaperPositionMode','auto');
eval(['print -depsc2 ',t1str,'_',t2str,'_dB_phi_poster',Cstr,'.eps']);