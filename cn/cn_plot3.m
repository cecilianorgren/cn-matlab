function [bB3 bB4] = cn_plot3(A,B,PA,PB)  

if size(A,1)>size(B,1); B=irf_resamp(B,A);
elseif size(A,1)<size(B,1); A=irf_resamp(A,B); 
end
nt=size(A,1);

%PA=[-5,5,-15];
%PB=-[-5,5,-15];

%dz=10;
%dz=1;dy=0;dx=-0.5;
%% Spacecraft positions

%t12=(t1+t2)/2;
%c_eval('gsePos?=cn_toepoch(t12,Pos?);',3:4);
% Transforming to field aligned system        
%c_Eval('Pos?=M*gsePos?(2:4)'';',3:4);
%Pos0=(Pos3+Pos4)/2;
%c_eval('Pos?c=Pos?-Pos0;',3:4);
%boxside=max([Pos3c(1) Pos3c(2) Pos4c(1) Pos4c(2)]);
%plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
%zdist=abs(Pos4c(3)-Pos3c(3));
%cn_mag(Pos3c-Pos4c);
xlabel('x_{B} [km]'); ylabel('y_{B} [km]');

figure('name','Fields');
irf_plot({A,B},'comp')
%% Plot Eperp and deltaB
% make position vector;
figure('name','Topology: B');
X=[repmat(PA(1),nt,1) repmat(PB(1),nt,1)]+repmat(dx*((1:nt)'),1,2);
Y=[repmat(PA(2),nt,1) repmat(PB(2),nt,1)]+repmat(dy*((1:nt)'),1,2);
Z=[repmat(PA(3),nt,1) repmat(PB(3),nt,1)]+repmat(dz*((1:nt)'),1,2);
%for k=1:nt
    %plot(Pos3c(1),Pos3c(2),'go',Pos4c(1),Pos4c(2),'bv'); hold on;
    %plot([Px3 Px4],[Py3 Py4],'color',line,'linestyle','--')
    quiver3(X(:,1),Y(:,1),Z(:,1),A(:,2),A(:,3),A(:,4),'color',[1 0 0])
    hold on;
    quiver3(X(:,2),Y(:,2),Z(:,2),B(:,2),B(:,3),B(:,4),'color',[0 1 0])
    %quiver3(Px3,Py3,bE3(l,2)*s,bE3(l,3)*s,0,'color',green)
    
    %Px3=Px3;%-dx;
    %Py3=Py3-dy;
    %Px4=Px4;%-dx;
    %Py4=Py4-dy;
    %pause(0.01)
%end
xlabel('x');ylabel('y');zlabel('z');

if 0
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
end
