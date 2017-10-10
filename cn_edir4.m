function [bB3 bB4] = cn_edir(t1,t2,ta,tb,E30,E40,ExB30,ExB40,gsmB3,gsmB4,Pos3,Pos4,vph,M,flh,n_hat,Te,t0,ACDC,scale)  
%% Find E AC and ExB DC
c_eval('E?=cn_toepoch(t1,t2,E?0);',3:4);
c_eval('ExB?=cn_toepoch(t1,t2,ExB?0);',3:4);
dt=E3(2,1)-E3(1,1);
fs=1/dt;
flow=flh*0.4;
fhigh=180;
% Taking a slightly larger interval to make the filtering
t11=t1;
t11(5)=t11(5)-1;
t22=t2;
t22(5)=t22(5)+1;

c_eval('E?2=cn_toepoch(t11,t22,E?0);',3:4)
c_eval('E?DC=irf_filt(E?2,0,flow,fs,3);',3:4);
c_eval('E?AC=[E?2(:,1) E?2(:,2:4)-E?DC(:,2:4)];',3:4)
c_eval('E?AC=cn_toepoch(ta,tb,E?AC);',3:4)
c_eval('bE?AC=M*E?AC(:,2:4)'';',3:4);
c_eval('bE?=[E?AC(:,1) bE?AC''];',3:4);

c_eval('bEt?=M*E?(:,2:4)'';',3:4);
c_eval('bEt?=[E?(:,1) bEt?''];',3:4);

c_eval('ExB?2=cn_toepoch(t11,t22,ExB?0);',3:4)
c_eval('ExB?DC=irf_filt(ExB?2,0,flow,fs,3);',3:4);
c_eval('ExB?DC=cn_toepoch(ta,tb,ExB?DC);',3:4)



switch ACDC
    case 'AC'
        c_eval('bE?=bE?;',3:4);
    case 'tot'        
        c_eval('bE?=bEt?;',3:4);
end

%ExB=irf_add(1,ExB3DC(:,2:4),1,ExB3DC(:,2:4))/2;
ExB=irf_add(0.5,ExB3DC,0.5,ExB4DC);
%ExB=sum(ExB,1)/size(ExB,1);
bExB=M*ExB(:,2:4)';
%% Find B AC 
c_eval('gseB?=irf_gse2gsm(gsmB?,-1);',3:4);
c_eval('B?=cn_toepoch(t11,t22,gseB?);',3:4);

c_eval('B?AC=irf_filt(B?,flow,0,fs,3);',3:4);
c_eval('B?AC=cn_toepoch(ta,tb,B?AC);',3:4)
c_eval('bB?AC=M*B?AC(:,2:4)'';',3:4);
c_eval('bB?=[B?AC(:,1) bB?AC''];',3:4);
%% Get E and B in same length
    bB3=irf_resamp(bB3,bE3);
    bB4=irf_resamp(bB4,bE3);
    bE4=irf_resamp(bE4,bE3);
%%
% size(bExB)
% size(ExB30)
bExB=[ExB(:,1) bExB'];
%figure;irf_plot(bExB);

%% Spacecraft positions
figure('name','Topology: E perp');
t12=(t1+t2)/2;
c_eval('gsePos?=cn_toepoch(t12,Pos?);',3:4);
% Transforming to field aligned system        
c_Eval('Pos?=M*gsePos?(2:4)'';',3:4);
Pos0=(Pos3+Pos4)/2;
c_eval('Pos?c=Pos?-Pos0;',3:4);
boxside=max([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]);
plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
zdist=abs(Pos4c(3)-Pos3c(3));
cn_mag(Pos3c-Pos4c);
xlabel('x_{B} [km]'); ylabel('y_{B} [km]');

%% Plot Eperp and deltaB
red=[1 0 0];
blue=[0 0.5 1];
green=[0 1 0.5];
line=[0.7 0.7 0];
max(E3(:,2))
s=scale/max(E3(:,2))
dy=-dt*vph;
%bE3=irf_resamp(bE3,bB3);
Px3=Pos3c(1);
Py3=Pos3c(2);
Px4=Pos4c(1);
Py4=Pos4c(2);
for l=1:1:size(E3AC,1) 
    dx=dt*bExB(l,2);
    if bB3(l,4) > 0
        plot(Px3,Py3,'ro');
    else
        plot(Px3,Py3,'rx');
    end
    if bB4(l,4) > 0
        plot(Px4,Py4,'ro');
    else
        plot(Px4,Py4,'rx');
    end
    plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
    plot([Px3 Px4],[Py3 Py4],'color',line,'linestyle',':')
    quiver(Px4,Py4,bE4(l,2)*s,bE4(l,3)*s,0,'color',blue)
    quiver(Px3,Py3,bE3(l,2)*s,bE3(l,3)*s,0,'color',green)
    
    Px3=Px3-dx;
    Py3=Py3-dy;
    Px4=Px4-dx;
    Py4=Py4-dy;
    %pause(0.01)
end
if 0 % Add length scale for the E-arrows
    onekm=1;
    text(10,10,'10 km corresponds to ')
end
if 0
for l=1:1:size(E3AC,1) 
    dx=-dt*bExB(l,3);
    if bB3(l,4) > 0
        plot(Pos3c(1)-(l-1)*dx,Pos3c(2)-(l-1)*dy,'ro');
    else
        plot(Pos3c(1)-(l-1)*dx,Pos3c(2)-(l-1)*dy,'rx');
    end
    if bB4(l,4) > 0
        plot(Pos4c(1)-(l-1)*dx,Pos4c(2)-(l-1)*dy,'ro');
    else
        plot(Pos4c(1)-(l-1)*dx,Pos4c(2)-(l-1)*dy,'rx');
    end
    %disp('dB3 = ',[num2str(dB3(l,2:4)),'  dB4 = ',num2str(dB4(l,2:4))])
    %E3AC(l,2:4);
    %c_eval('bE?=M*E?AC(l,2:4)'';',3:4)
    plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
    plot([Pos3c(1)-(l-1)*dx Pos4c(1)-(l-1)*dx],...
        [Pos3c(2)-(l-1)*dy Pos4c(2)-(l-1)*dy],'color',line,'linestyle',':')
    quiver(Pos4c(1)-(l-1)*dx,Pos4c(2)-(l-1)*dy,bE4(l,2)*s,bE4(l,3)*s,0,'color',blue)
    quiver(Pos3c(1)-(l-1)*dx,Pos3c(2)-(l-1)*dy,bE3(l,2)*s,bE3(l,3)*s,0,'color',green)
    pause(0.01)
end
end
if 0 % Add timescale to the right of plot
    yticks=get(hcf,'ytick')
end
%% Adding info to figure
time1=[num2str(t1(1)),'-',num2str(t1(2)),'-',num2str(t1(3)),...
    ' ',num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)),...
    '.',num2str(t1(7))];
time2=[num2str(t2(1)),'-',num2str(t2(2)),'-',num2str(t2(3)),...
    ' ',num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)),...
    '.',num2str(t2(7))];

title([time1,' to ',time2,'    v=',num2str(vph,'%.0f'),' km/s']);
legend('C3','C4');
set(gcf,'PaperPositionMode','auto'); 
ax1=gca;
%% Plot bE timeseries
if 0
    figure('name','bE timeseries');
    h=irf_plot(5);
    isub=1;

    if 1 % Plot bE 
        c_eval('bE?shift=bE?;',3:4);
        c_eval('bE?shift(:,1)=bE?shift(:,1)-t0;',3:4);
        hca=h(isub);isub=isub+1;
        irf_plot(hca,bE3shift(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,bE4(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{\perp k}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        irf_zoom(hca,'y');  
        irf_zoom(hca,'x',[bE3(1,1) bE3(end,1)])   

        hca=h(isub);isub=isub+1;
        irf_plot(hca,bE3shift(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,bE4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{|| k}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);   
        irf_zoom(hca,'y'); 
        irf_zoom(hca,'x',[bE3(1,1) bE3(end,1)])
    end
    if 1 % plot delta B C3 C4    
        hca=h(isub);isub=isub+1;
        irf_plot(hca,bB3); hold(hca,'on');
        ylabel(hca,'C3    \delta B [nT]');
        irf_zoom(hca,'y'); 
        irf_zoom(hca,'x',[bE3(1,1) bE3(end,1)])
        irf_legend(hca,{'x','y','z'},[0.02 0.05]); 

        hca=h(isub);isub=isub+1;
        irf_plot(hca,bB4); hold(hca,'on');
        ylabel(hca,'C4    \delta B [nT]');
        irf_zoom(hca,'y'); 
        irf_zoom(hca,'x',[bE4(1,1) bE4(end,1)])
        irf_legend(hca,{'x','y','z'},[0.02 0.05]); 
    end
    if 1 % plot delta B z    
        hca=h(isub);isub=isub+1;
        irf_plot(hca,bB3(:,[1 4]),'g'); hold(hca,'on');
        irf_plot(hca,bB4(:,[1 4]),'b'); hold(hca,'on');
        ylabel(hca,'C3    \delta B_Z [nT]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        irf_zoom(hca,'y'); 
        irf_zoom(hca,'x',[bE3(1,1) bE3(end,1)])
    end

    %axis equal
    if 0
    annotation('textbox',...
    [0.704571428571426 0.183333333333336 0.241857142857146 0.345238095238099],...
    'String',{['\rho_e = ',num2str(r_e,'%.1f'),'km'],...
        ['f_{lh} = ',num2str(flh,'%.1f'),'Hz'],...
        ['d_{34,||} = ',num2str(zdist,'%.1f'),'km']},...           
    'FitBoxToText','off','LineStyle','none',...
    'BackgroundColor',[1 1 1]);
    end
        set(gcf,'PaperPositionMode','auto'); 
end