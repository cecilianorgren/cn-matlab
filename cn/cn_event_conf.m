function f=cn_event_conf(t1,t2,B3_raw,B4_raw,Vi3_raw,Vi4_raw,ExB3_raw,ExB4_raw,Pos3_raw,Pos4_raw,Te,Ti,Ne,Ni,Ne3,PeaNe3,n_hat,t1str,t2str,way)

switch way
    case 'exjobb'
tstart=t1;
tstop=t2;
figure('name','Event configuration')
%% Finding average bhat under the time interval
c_eval('B?=cn_toepoch(t1,t2,B?_raw);',3:4);
[b_av b_hat b_mag]=cn_hat(B3,B4);
%% Defining coordinate system
z=b_hat;
x=cn_cross(cn_cross(z,n_hat),z); % close to BL-normal
y=cn_cross(z,x); % close to drift direction within boundary layer
M=[x;y;z]; % transformation matrix
%% Taking average Vi under the time interval
c_eval('Vi?=cn_toepoch(t1,t2,Vi?_raw);',3:4);
[vi_av vi_hat vi_mag]=cn_hat(Vi3,Vi4);
%% Taking average ExB-velocity under the time interval
c_eval('ExB?=cn_toepoch(t1,t2,ExB?_raw);',3:4);
[ExB_av ExB_hat ExB_mag]=cn_hat(ExB3,ExB4);
%% Spacecraft positions
t12=fix((t1+t2)/2);
c_eval('gsePos?=cn_toepoch(t12,Pos?_raw);',3:4);

% Transforming to field aligned system        
Pos3=M*gsePos3(2:4)';
Pos4=M*gsePos4(2:4)';
Pos0=(Pos3+Pos4)/2;
c_eval('Pos?c=Pos?-Pos0;',3:4);
%axes('Position',[-7 -7 7*2 7*2])

boxside=max([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]);
plot(Pos3c(1),Pos3c(2),'go'); hold on;
plot(Pos4c(1),Pos4c(2),'bv'); hold on;

normdist=abs(Pos4c(1)-Pos3c(1));
zdist=abs(Pos4c(3)-Pos3c(3));
kdist=abs(Pos4c(2)-Pos3c(2));
cn_mag(Pos3c-Pos4c);

xlabel('x_{B} [km]'); ylabel('y_{B} [km]');       
axis([-boxside boxside -boxside boxside]*1.1)
axis square       
%% Plot Viav
mvi=(M*vi_av')';
mvi_hat=(M*vi_hat')';
V_par=mvi(3);
V_per=sqrt(mvi(1)^2+mvi(2)^2);
quiver(0,0,mvi_hat(1),mvi_hat(2),4,'color',[0.8 0.8 0]); hold on; 
%% Plot ExBav
%ExB=M*ExBav(2:4)';
%ExB_per=sqrt(ExB(1)^2+ExB(2)^2);
%ExB_par=ExB(3);
%ExBhat=ExB/cn_mag(ExB);
mExB_hat=(M*ExB_hat')';
mExB_av=(M*ExB_av')';
mExB_per=sqrt(mExB_av(1)^2+mExB_av(2)^2);
mExB_par=mExB_av(3);
quiver(0,0,mExB_hat(1),mExB_hat(2),4,'color',[1 0.8 0]); hold on; 

if 0 % plot xyz gse
%% Plot "towards earth" and BL-normal
xb=M*[1 0 0]';
quiver(0,0,xb(1),xb(2),1,'b'); hold on;
yb=M*[0 1 0]';
quiver(0,0,yb(1),yb(2),1,'c'); hold on;
zb=M*[0 0 1]';
quiver(0,0,zb(1),zb(2),1,'r'); hold on;
end 
%% Plot ISR2x and ISR2y
isr2x=c_coord_trans('DSI','GSE',[cn_toepoch(fix((t1+t2)/2)) 1 0 0],'CL_ID',3);              
isr2y=c_coord_trans('DSI','GSE',[cn_toepoch(fix((t1+t2)/2)) 0 1 0],'CL_ID',3);
isr2xb=M*isr2x(2:4)';
isr2yb=M*isr2y(2:4)';
quiver(0,0,isr2xb(1),isr2xb(2),2,'b'); hold on;
quiver(0,0,isr2yb(1),isr2yb(2),2,'g'); hold on;
%% Boundary normal
bn_hat=M*n_hat';
quiver(0,0,1,0,4); hold on;       
%% Calculating ro_e, ro_i, L_n, flh, veth, vith (etc?) 
Tez=cn_toepoch(t1,Te);
Te_av=cn_average2(Tez,1);
Tiz=cn_toepoch(t1,t2,Ti);
Ti_av=cn_average2(Tiz,1);

% Electron gyroradius
r_e=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Roe'); % in m
r_e=r_e/1000; % in km      

% Ion gyroradius 
r_i=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Rop'); % in m
r_i=r_i/1000; % km

% Lower hybrid frequency
flh=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Flh');

% Ion thermal velocity     
Vith=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Vtp'); % m/s        
Vith=Vith/1000; % km/s

% Electron thermal velcoity
Veth=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Vte'); % m/s        
Veth=Veth/1000; % km/s

% Gradient scale length
Ln=cn_ln(t1,t2,ExB3_raw,flh,Ne3,PeaNe3,M);

% Adding Debye length
%Nem3=Ne*1.0000e-06;
%Deb=7.430*sqrt(Teav/Nem3);
%Deb2=sqrt(200*Teav/Ne);
%% Transforming kvec to gse for returning
vd_hat=(M*n_hat')';        
%% B-to-spin-plane angle
%cosalfa=cn_scalar(Bhat,c_coord_trans('DSI','GSE',[cn_toepoch((t1+t2)/2) [1 1 0]/sqrt(2)],'CL_ID',3));
%alfa=acosd(cosalfa(2));
%% Adding title, legends and info box to figure
time1=[num2str(t1(1)),'-',num2str(t1(2)),'-',num2str(t1(3)),...
    ' ',num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)),...
    '.',num2str(t1(7))];
time2=[num2str(t2(1)),'-',num2str(t2(2)),'-',num2str(t2(3)),...
    ' ',num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)),...
    '.',num2str(t2(7))];

title([time1,' to ',time2]);%,...
    %'\newline \rho_e=',num2str(r_e),'km   v_{||}=',num2str(V_par),...
    %'km/s   v_{\perp}=',num2str(V_per),'km/s   d_{34,||}=',num2str(zdist),'km',...
    %'\newline \lambda_D=',num2str(Deb),'km'])
legend('C3','C4','v_{\perp,i}',...
    'v_{\perp,ExB}','x^{ISR2}','y^{ISR2}','BL-normal',...
    'Location','NorthEastOutside')


annotation('textbox',...
[0.734928571428569 0.247619047619051 0.241857142857146 0.345238095238099],...
'String',{['\rho_e = ',num2str(r_e,'%.1f'),' km'],...
    ['f_{lh} = ',num2str(flh,'%.1f'),' Hz'],...
    ['\Delta z_{C3,C4} = ',num2str(zdist,'%.0f'),' km'],...
    ['v_{||} = ',num2str(V_par,'%.0f'),' km/s'],...
    ['v_{\perp} = ',num2str(V_per,'%.0f'),' km/s'],...
    ['v_{ExB,\perp} = ',num2str(mExB_per,'%.0f'),' km/s'],...
    ['v_{th,i} = ',num2str(Vith,'%.0f'),' km/s'],... 
    ['L_n/ \rho_{ / i} = ',num2str(Ln/r_i,'%.1f'),'']... 
    %,['v_{th,e} = ',num2str(Veth,'%.1f'),'km/s']...
    },...           
'FitBoxToText','off','LineStyle','none',...
'BackgroundColor',[1 1 1]);
%% Returning values
mb_av=(M*b_av')';
f= {b_hat M mvi vi_hat mb_av [0 1 0] flh Vith Ti_av Te_av normdist kdist r_i zdist r_e Ln}; 


set(gcf,'PaperPositionMode','auto'); 
eval(['print -depsc2 ',t1str,'_',t2str,'_plot.eps']);
    case 'poster'
        
tstart=t1;
tstop=t2;
figure('name','Event configuration')
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');

%% Finding average bhat under the time interval
c_eval('B?=cn_toepoch(t1,t2,B?_raw);',3:4);
[b_av b_hat b_mag]=cn_hat(B3,B4);
%% Defining coordinate system
z=b_hat;
x=cn_cross(cn_cross(z,n_hat),z); % close to BL-normal
x=x/cn_mag(x);
y=cn_cross(z,x); % close to drift direction within boundary layer
y=y/cn_mag(y);
M=[x;y;z]; % transformation matrix
%% Taking average Vi under the time interval
c_eval('Vi?=cn_toepoch(t1,t2,Vi?_raw);',3:4);
[vi_av vi_hat vi_mag]=cn_hat(Vi3,Vi4);
%% Taking average ExB-velocity under the time interval
c_eval('ExB?=cn_toepoch(t1,t2,ExB?_raw);',3:4);
[ExB_av ExB_hat ExB_mag]=cn_hat(ExB3,ExB4);
%% Calculating ro_e, ro_i, L_n, flh, veth, vith (etc?) 
Tez=cn_toepoch(t1,Te);
Te_av=cn_average2(Tez,1);
Tiz=cn_toepoch(t1,t2,Ti);
Ti_av=cn_average2(Tiz,1);

% Electron gyroradius
r_e=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Roe'); % in m
r_e=r_e/1000; % in km      

% Ion gyroradius 
r_i=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Rop'); % in m
r_i=r_i/1000; % km

% Lower hybrid frequency
flh=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Flh');

% Ion thermal velocity     
Vith=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Vtp'); % m/s        
Vith=Vith/1000; % km/s

% Electron thermal velocity
Veth=irf_plasma_calc(b_mag,Ni,0,Te_av,Ti_av,'Vte'); % m/s        
Veth=Veth/1000; % km/s

% Gradient scale length
Ln=cn_ln(t1,t2,ExB3_raw,flh,Ne3,PeaNe3,M);


%% Spacecraft positions
t12=fix((t1+t2)/2);
c_eval('gsePos?=cn_toepoch(t12,Pos?_raw);',3:4);

% Transforming to field aligned system        
Pos3=M*gsePos3(2:4)';
Pos4=M*gsePos4(2:4)';
Pos0=(Pos3+Pos4)/2;
c_eval('Pos?c=Pos?-Pos0;',3:4);
%axes('Position',[-7 -7 7*2 7*2])

boxside=max(abs([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]));
plot(Pos3c(1),Pos3c(2),'go','markersize',13,'linewidth',1.5); hold on;
plot(Pos4c(1),Pos4c(2),'bv','markersize',13,'linewidth',1.5); hold on;
c3=text(Pos3c(1),Pos3c(2),'    C3','color','k');
c4=text(Pos4c(1),Pos4c(2),'    C4','color','k');

normdist=abs(Pos4c(1)-Pos3c(1));
zdist=Pos4c(3)-Pos3c(3);
kdist=abs(Pos4c(2)-Pos3c(2));
cn_mag(Pos3c-Pos4c);

xlabel('x   [km]','fontsize',16); ylabel('y   [km]','fontsize',16);       
axis([-boxside boxside -boxside boxside]*1.1)
axis square       

set(gcf,'PaperPositionMode','auto'); 

%% Boundary normal
bn_hat=M*n_hat';
ncol=[1 0 0.7];
q1=quiver(0,0,1,0,4,'color',ncol); hold on;  

nstr(1)={'boundary'};
nstr(2)={'normal'};
text(3.3,-0.75,nstr,'color',ncol);


%% Add B_0
bcol=[0.6 0 0];
bdot  = plot(-1,-1,'color',bcol,'linestyle','o','markersize',2,'linewidth',1);
bring = plot(-1,-1,'color',bcol,'linestyle','o','markersize',20,'linewidth',1.1);
btext = text(-0.5,-1,'  B','color',bcol);


%% Plot Viav
mvi=(M*vi_av')';
mvi_hat=(M*vi_hat')';
V_par=mvi(3);
V_per=sqrt(mvi(1)^2+mvi(2)^2);
q2=quiver(0,0,mvi_hat(1),mvi_hat(2),2*V_per/Vith,'color',[0.5 0.8 0.1]); hold on; 
vistr(1)={'ion velocity '};
vistr(2)={'from CIS'};
text(-1.2,3.7,vistr,'color',[0.5 0.8 0.1]);
%% Plot ExBav
mExB_hat=(M*ExB_hat')';
mExB_av=(M*ExB_av')';
mExB_per=sqrt(mExB_av(1)^2+mExB_av(2)^2);
mExB_par=mExB_av(3);
q3=quiver(0,0,mExB_hat(1),mExB_hat(2),2*mExB_per/Vith,'color',[1 0.6 0],'linewidth',1.1); hold on; 

text(2,3,'ExB-velocity','color',[1 0.6 0]);

%% Propagation direction
bn_hat=M*n_hat';
kcol=[0.5 0 1];
q4=quiver(0,0,0,975,2/Vith,'color',[0.5 0 1]); hold on;
%quiver(0,0,0,-1,4,'color',[0 0 1]); hold on;
kstr(1)={'wave'};
kstr(2)={'velocity'};
text(-2,2,kstr,'color',kcol);

%% Transforming kvec to gse for returning
vd_hat=(M*n_hat')';        
%% B-to-spin-plane angle
%cosalfa=cn_scalar(Bhat,c_coord_trans('DSI','GSE',[cn_toepoch((t1+t2)/2) [1 1 0]/sqrt(2)],'CL_ID',3));
%alfa=acosd(cosalfa(2));
%% Adding title, legends and info box to figure
time1=[num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6))];
time2=[num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6))];

title(['August 31 2007  ' time1,' to ',time2],'fontsize',16);%,...
    %'\newline \rho_e=',num2str(r_e),'km   v_{||}=',num2str(V_par),...
    %'km/s   v_{\perp}=',num2str(V_per),'km/s   d_{34,||}=',num2str(zdist),'km',...
    %'\newline \lambda_D=',num2str(Deb),'km'])
%legend('C3','C4','v_{\perp,i}',...
%    'v_{\perp,ExB}','x^{ISR2}','y^{ISR2}','BL-normal',...
%    'Location','NorthEastOutside')


infocolor=[0 0 0];
if 0 % black plot
    set(gcf,'defaultTextColor',[1 1 1]);
    
    infocolor=[1 1 1];
    bgcolor=[0 0 0];
    bgfcolor=[0 0 0];
    frame=[1 1 1];
    bgfcolor=[0 0 0];
    bgfcolor=[0 0 0];
    
    ylab=get(gca,'ylabel');
    set(ylab,'color',infocolor);
    xlab=get(gca,'xlabel');
    set(xlab,'color',infocolor);
    set(gcf,'color',bgfcolor)
    set(gca,'color',bgcolor)
    set(gca,'xcolor',frame)
    set(gca,'ycolor',frame)
    set(gcf, 'InvertHardCopy', 'off');
    
    set(c3,'color',[0 1 0]);
    set(c4,'color',[0 0 1]);
    set(bdot,'color',[1 0 0]);
    set(bring,'color',[1 0 0]);
    set(btext,'color',[1 0 0]);
end    
    set(q1,'LineWidth',1.1);
    set(q2,'LineWidth',1.1);
    set(q3,'LineWidth',1.1);
    set(q4,'LineWidth',1.1);


set(gcf,'defaultTextFontSize',14);

infostr(1)={['\rho_e = ',num2str(r_e,'%.1f'),' km']};
infostr(2)={['f_{lh} = ',num2str(flh,'%.1f'),' Hz']};
infostr(3)={['\Delta z_{C3,C4} = ',num2str(zdist,'%.0f'),' km']};
infostr(4)={['v_{||,CIS} = ',num2str(V_par,'%.0f'),' km/s']};
infostr(5)={['v_{\perp,CIS} = ',num2str(V_per,'%.0f'),' km/s']};
infostr(6)={['v_{ExB} = ',num2str(mExB_per,'%.0f'),' km/s']};
infostr(7)={['v_{th,i} = ',num2str(Vith,'%.0f'),' km/s']};
infostr(8)={['v_{wave} = 978 km/s']};
text(-boxside,-boxside,infostr,'color',infocolor,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');


% ann=annotation('textbox',...
% [0.23135714285714 0.140476190476195 0.207928571428574 0.305238095238099],...
% 'BackgroundColor',[0 0 1],'Color',[0 1 1],...
% 'String',{['\rho_e = ',num2str(r_e,'%.1f'),' km'],...
%     ['f_{lh} = ',num2str(flh,'%.1f'),' Hz'],...
%     ['\Delta z_{C3,C4} = ',num2str(zdist,'%.0f'),' km'],...
%     ['v_{||,CIS} = ',num2str(V_par,'%.0f'),' km/s'],...
%     ['v_{\perp,CIS} = ',num2str(V_per,'%.0f'),' km/s'],... %['v_{ExB,||} = ',num2str(mExB_par,'%.0f'),' km/s'],...
%     ['v_{ExB} = ',num2str(mExB_per,'%.0f'),' km/s'],...
%     ['v_{th,i} = ',num2str(Vith,'%.0f'),' km/s'],... 
%     ['v_{wave} = 978 km/s']...
%     %['L_n/ \rho_{i} = ',num2str(Ln/r_i,'%.1f'),'']... 
%     %,['v_{th,e} = ',num2str(Veth,'%.1f'),'km/s']...
%     },...           
% 'FitBoxToText','off','LineStyle','none',...
% 'BackgroundColor',[1 1 1]);
% get(ann)
%% Returning values
mb_av=(M*b_av')';
f= {b_hat M mvi vi_hat mb_av [0 1 0] flh Vith Ti_av Te_av normdist kdist r_i zdist r_e Ln}; 



eval(['print -depsc2 ',t1str,'_',t2str,'_plot_poster.eps']);
end