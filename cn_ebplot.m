function [bB3 bB4] = cn_edir(t1,t2,E30,E40,gsmB3,gsmB4,Pos3,Pos4,vph,M,flh,r_e,n_hat,Te,t0,ACDC,scale)  
%% Find E AC and ExB DC
x=M(1,:);
y=M(2,:);
z=M(3,:);

%M(1,:)=y;
%M(2,:)=x;
c_eval('E?=cn_toepoch(t1,t2,E?0);',3:4);
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
c_eval('E?AC=cn_toepoch(t1,t2,E?AC);',3:4)
c_eval('bE?AC=M*E?AC(:,2:4)'';',3:4);
c_eval('bE?=[E?AC(:,1) bE?AC''];',3:4);

c_eval('bEt?=M*E?(:,2:4)'';',3:4);
c_eval('bEt?=[E?(:,1) bEt?''];',3:4);


switch ACDC
    case 'AC'
        c_eval('bE?=bE?;',3:4);
    case 'tot'        
        c_eval('bE?=bEt?;',3:4);
end


%% Find B AC 
c_eval('gseB?=irf_gse2gsm(gsmB?,-1);',3:4);
c_eval('B?=cn_toepoch(t11,t22,gseB?);',3:4);

c_eval('B?AC=irf_filt(B?,flow,0,fs,3);',3:4);
c_eval('B?AC=cn_toepoch(t1,t2,B?AC);',3:4)
c_eval('bB?AC=M*B?AC(:,2:4)'';',3:4);
c_eval('bB?=[B?AC(:,1) bB?AC''];',3:4);
%% Get E and B in same length
    bB3=irf_resamp(bB3,bE3);
    bB4=irf_resamp(bB4,bE3);
    bE4=irf_resamp(bE4,bE3);

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
line=[0.7 0.7 0.7];
max(E3(:,2))
s=scale/max(E3(:,2))
dy=-dt*vph;
%bE3=irf_resamp(bE3,bB3);
Px3=Pos3c(1);
Py3=Pos3c(2);
Px4=Pos4c(1);
Py4=Pos4c(2);
for l=1:1:size(E3AC,1) 
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
    plot([Px3 Px4],[Py3 Py4],'color',line,'linestyle','--')
    quiver(Px4,Py4,bE4(l,2)*s,bE4(l,3)*s,0,'color',blue)
    quiver(Px3,Py3,bE3(l,2)*s,bE3(l,3)*s,0,'color',green)
    
    Px3=Px3;%-dx;
    Py3=Py3-dy;
    Px4=Px4;%-dx;
    Py4=Py4-dy;
    %pause(0.01)
end

%% Adding info to figure
time1=[num2str(t1(1)),'-',num2str(t1(2)),'-',num2str(t1(3)),...
    ' ',num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)),...
    '.',num2str(t1(7))];
time2=[num2str(t2(1)),'-',num2str(t2(2)),'-',num2str(t2(3)),...
    ' ',num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)),...
    '.',num2str(t2(7))];

title(['Wave perpendicular electric and parallel magnetic field']);
legend('C3','C4');
set(gcf,'PaperPositionMode','auto'); 
ax1=gca;
