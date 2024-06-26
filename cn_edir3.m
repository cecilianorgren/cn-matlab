function cn_edir(t1,t2,E3_raw,E4_raw,B3_raw,B4_raw,Vi3_raw,Vi4_raw,ExB3_raw,ExB4_raw,Pos3_raw,Pos4_raw,Te,Ti,Ne,Ni,n_hat,t1str,t2str,whichplot)

%%% Taking average B under the time interval
c_eval('B?=cn_toepoch(t1,t2,B?_raw);',3:4);
c_eval('B?av=cn_average(B?);',3:4);
Bav=cn_average(B3av,B4av);
magBav=cn_mag(Bav);
Bhat=[Bav(1) Bav(2:4)/magBav(2)]; % Bhat     

%% Defining coordinate system
z=Bhat(2:4);
x=cn_cross(cn_cross(z,n_hat),z); % close to BL-normal
y=cn_cross(z,x); % close to drift direction within boundary layer

M=[x;y;z];         

%% Spacecraft positions
t12=(t1+t2)/2;
c_eval('gsePos?=cn_toepoch(t12,Pos?_raw);',3:4);

% Transforming to field aligned system        
Pos3=M*gsePos3(2:4)';
Pos4=M*gsePos4(2:4)';
Pos0=(Pos3+Pos4)/2;
c_eval('Pos?c=Pos?-Pos0;',3:4);

boxside=max([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]);

plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;

zdist=abs(Pos4c(3)-Pos3c(3));
cn_mag(Pos3c-Pos4c);
xlabel('x_{B} [km]'); ylabel('y_{B} [km]');
%axis([-boxside boxside -boxside*5 boxside]*10)
%axis equal
%axis square

if 0 % Transform E
step=2;
c_eval('bE?=zeros(ceil(size(E3,1)/step),4);',3:4);
for l=1:step:length(E3); bE3(ceil(l/step),2:4)=M*E3(l,2:4)'; end
for l=1:step:length(E4); bE4(ceil(l/step),2:4)=M*E4(l,2:4)'; end
end
%% Plot E-direction
c_eval('E?=cn_toepoch(t1,t2,E?_raw);',3:4);
s=50/max(E3(:,2));
dr=0.5
for l=1:5:size(E3)  
    E3(l,2:4);
    c_eval('bE?=M*E?(l,2:4)'';',3:4)
    plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
    quiver(Pos4c(1),Pos4c(2)-l*dr,bE4(1)*s,bE4(2)*s,0)
    quiver(Pos3c(1),Pos3c(2)-l*dr,bE3(1)*s,bE3(2)*s,0)
    pause(0.1)
end
if 0% Plot ExBav
ExB=M*ExBav(2:4)';
ExB_per=sqrt(ExB(1)^2+ExB(2)^2);
ExB_par=ExB(3);
ExBhat=ExB/cn_mag(ExB);
quiver(0,0,ExBhat(1),ExBhat(2),4,'color',[1 0.8 0]); hold on; 
end
if 0 % plot xyz gse
    
%% Plot "towards earth" and BL-normal
xb=M*[1 0 0]';
quiver(0,0,xb(1),xb(2),1,'b'); hold on;
yb=M*[0 1 0]';
quiver(0,0,yb(1),yb(2),1,'c'); hold on;
zb=M*[0 0 1]';
quiver(0,0,zb(1),zb(2),1,'r'); hold on;
end 
if 0 % Plot ISR2x and ISR2y
isr2x=c_coord_trans('DSI','GSE',[cn_toepoch((t1+t2)/2) 1 0 0],'CL_ID',3);              
isr2y=c_coord_trans('DSI','GSE',[cn_toepoch((t1+t2)/2) 0 1 0],'CL_ID',3);
isr2xb=M*isr2x(2:4)';
isr2yb=M*isr2y(2:4)';
quiver(0,0,isr2xb(1),isr2xb(2),2,'b'); hold on;
quiver(0,0,isr2yb(1),isr2yb(2),2,'g'); hold on;

%ExB
nhatb=M*n_hat';
quiver(0,0,nhatb(1),nhatb(2),4); hold on;       
end

%% Adding electron gyro radius
Tezoom=cn_toepoch(t1,Te);
Teav=cn_average(Tezoom);
r_e=irf_plasma_calc(Bav(5),0.1,0,Teav,2500,'Roe'); % in m
r_e=r_e/1000; % in km      

% Lower hybrid frequency
flh=irf_plasma_calc(cn_mag(Bav(2:4)),Ni,0,Teav,3000,'Flh');


%% B-to-spin-plane angle
cosalfa=cn_scalar(Bhat,c_coord_trans('DSI','GSE',[cn_toepoch((t1+t2)/2) [1 1 0]/sqrt(2)],'CL_ID',3));
alfa=acosd(cosalfa(2))

%% Adding info to figure
time1=[num2str(t1(1)),'-',num2str(t1(2)),'-',num2str(t1(3)),...
    ' ',num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)),...
    '.',num2str(t1(7))];
time2=[num2str(t2(1)),'-',num2str(t2(2)),'-',num2str(t2(3)),...
    ' ',num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)),...
    '.',num2str(t2(7))];

title([time1,' to ',time2,'  \alpha =',num2str(alfa)]);
legend('C3','C4',...
    'BL-normal',...
    'Location','NorthEastOutside');


annotation('textbox',...
[0.704571428571426 0.183333333333336 0.241857142857146 0.345238095238099],...
'String',{['\rho_e = ',num2str(r_e,'%.1f'),'km'],...
    ['f_{lh} = ',num2str(flh,'%.1f'),'Hz'],...
    ['d_{34,||} = ',num2str(zdist,'%.1f'),'km']},...           
'FitBoxToText','off','LineStyle','none',...
'BackgroundColor',[1 1 1]);
 
   