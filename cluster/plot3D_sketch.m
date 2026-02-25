c_eval('gsex?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 1 0 0],''CL_ID'',?);',3:4);
c_eval('gsey?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 0 1) 0],''CL_ID'',?);',3:4);

c_eval('gseEx?=c_coord_trans(''dsi'',''gse'',[Ee3(1:2) 0 0],''CL_ID'',?);',3:4);
c_eval('gseEy?=c_coord_trans(''dsi'',''gse'',[Ee3(1) 0 Ee3(3) 0],''CL_ID'',?);',3:4);

%%
figure
s=30;
m=-110;
plot3(s*[0 1],s*[0 0],s*[0 0],'--c'); hold on;
plot3(s*[0 0],s*[0 1],s*[0 0],'--c'); hold on;
plot3(s*[0 0],s*[0 0],s*[0 1],'--c'); hold on;


plot3([0 0],[0 11.6],[0 -0.3],'b'); hold on;
plot3([0 52.8],[0 -0.2],[0 -6.8],'g'); hold on;
plot3(m*[0 -0.8],m*[0 -0.1],m*[0 -0.6],'r'); hold on;

%%
cn_scalar(Bh,gsex3(2:4))
cn_mag(gseEx3(2:4))

%%
ee3=cn_toepoch(t1,t2,diE3);
ee3n=bsxfun(@rdivide,mean(ee3(:,2:4),1),sqrt(sum(mean(ee3(:,2:4),1).^2,2)));
ee3p=ee3n;ee3p(:,4)=0;
e3n=irf_norm(ee3);

ee4=cn_toepoch(t1,t2,diE4);
ee4n=bsxfun(@rdivide,mean(ee4(:,2:4),1),sqrt(sum(mean(ee4(:,2:4),1).^2,2)));
ee4p=ee4n;ee4p(:,4)=0;
e4n=irf_norm(ee4);

bb3=cn_toepoch(t1,t2,diB3);
bb3n=bsxfun(@rdivide,mean(bb3(:,2:4),1),sqrt(sum(mean(bb3(:,2:4),1).^2,2)));
b3n=irf_norm(bb3);

bb4=cn_toepoch(t1,t2,diB4);
bb4n=bsxfun(@rdivide,mean(bb4(:,2:4),1),sqrt(sum(mean(bb4(:,2:4),1).^2,2)));
b4n=irf_norm(bb4);
%%
figure;
h=plot3(0,0,0,0,0,0);
ksum=0;
for k=1:100:size(bb4n,1)
    ksum=ksum+k;
    cn_plot3d(h,0.5*[dx(k,2) dy(k,2) ksum/100],b3n(k,2:4),'b');
    cn_plot3d(h,0.5*[-dx(k,2) -dy(k,2) ksum/100],b4n(k,2:4),'b');
end
%%
figure(19);
set(gcf,'defaultTextFontSize',14);
s=1;
% DSI
plot3(s*[0 1],s*[0 0],s*[0 0],'-c'); hold on;
plot3(s*[0 0],s*[0 1],s*[0 0],'-c'); hold on;
plot3(s*[0 0],s*[0 0],s*[0 1],'-c'); hold on;
text(1,0,0,'X_{DSI}','color',[0 1 1]);
text(0,1,0,'Y_{DSI}','color',[0 1 1]);
text(0,0,1,'Z_{DSI}','color',[0 1 1])

% B
plot3(s*[0 bb3n(1,1)],s*[0 bb3n(1,2)],s*[0 bb3n(1,3)],'-b'); hold on;
plot3(s*[0 bb4n(1,1)],s*[0 bb4n(1,2)],s*[0 bb4n(1,3)],'-r'); hold on;
text(bb3n(1),bb3n(2),bb3n(3),'BC3','color',[0 0 1])
text(bb4n(1),bb4n(2),bb4n(3),'BC4','color',[1 0 0])

% E
plot3(s*[0 ee3n(1,1)],s*[0 ee3n(1,2)],s*[0 ee3n(1,3)],'-b'); hold on;
plot3(s*[0 ee4n(1,1)],s*[0 ee4n(1,2)],s*[0 ee4n(1,3)],'-r'); hold on;
text(ee3n(1),ee3n(2),ee3n(3),'EC3','color',[0 0 1])
text(ee4n(1),ee4n(2),ee4n(3),'EC4','color',[1 0 0])
% E projection
plot3(s*[0 ee3n(1,1)],s*[0 ee3n(1,2)],s*[0 0],'--b'); hold on;
plot3(s*[0 ee4n(1,1)],s*[0 ee4n(1,2)],s*[0 0],'--r'); hold on;
text(ee3n(1),ee3n(2),0,'EC3','color',[0 0 1])
text(ee4n(1),ee4n(2),0,'EC4','color',[1 0 0])

% GSM
xgsm=[bb3(1,1) 1 0 0];
ygsm=[bb3(1,1) 0 1 0];
zgsm=[bb3(1,1) 0 0 1];
dixgsm=c_coord_trans('gsm','dsi',xgsm,'cl_id',4);
diygsm=c_coord_trans('gsm','dsi',ygsm,'cl_id',4);
dizgsm=c_coord_trans('gsm','dsi',zgsm,'cl_id',4);
plot3(s*[0 dixgsm(1,2)],s*[0 dixgsm(1,3)],s*[0 dixgsm(1,4)],'-g'); hold on;
plot3(s*[0 diygsm(1,2)],s*[0 diygsm(1,3)],s*[0 diygsm(1,4)],'-g'); hold on;
plot3(s*[0 dizgsm(1,2)],s*[0 dizgsm(1,3)],s*[0 dizgsm(1,4)],'-g'); hold on;
text(dixgsm(1,2),dixgsm(1,3),dixgsm(1,4),'X_{GSM}','color',[0 1 0])
text(diygsm(1,2),diygsm(1,3),diygsm(1,4),'Y_{GSM}','color',[0 1 0])
text(dizgsm(1,2),dizgsm(1,3),dizgsm(1,4),'Z_{GSM}','color',[0 1 0])

% FAC
aaa = [bb3(1,1) -0.81 0.54 0.25];
bbb = [bb3(1,1) -0.09 -0.52 0.85];
ccc = [bb3(1,1) 0.59 0.66 0.46];
aaagsm=c_coord_trans('gsm','dsi',aaa,'cl_id',4);
bbbgsm=c_coord_trans('gsm','dsi',bbb,'cl_id',4);
cccgsm=c_coord_trans('gsm','dsi',ccc,'cl_id',4);
plot3(s*[0 aaagsm(1,2)],s*[0 aaagsm(1,3)],s*[0 aaagsm(1,4)],'-k'); hold on;
plot3(s*[0 bbbgsm(1,2)],s*[0 bbbgsm(1,3)],s*[0 bbbgsm(1,4)],'-k'); hold on;
plot3(s*[0 cccgsm(1,2)],s*[0 cccgsm(1,3)],s*[0 cccgsm(1,4)],'-k'); hold on;
text(aaagsm(1,2),aaagsm(1,3),aaagsm(1,4),'X_{FAC}','color',[0 0 0])
text(bbbgsm(1,2),bbbgsm(1,3),bbbgsm(1,4),'Y_{FAC}','color',[0 0 0])
text(cccgsm(1,2),cccgsm(1,3),cccgsm(1,4),'Z_{FAC}','color',[0 0 0])


axis equal
set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
hold(gca,'off')
%% Angle of yfac above spin plane
x=aaagsm(2:4);
xang=atand(abs(x(3))/sqrt(x(1)^2+x(2)^2));

y=bbbgsm(2:4);
yang=atand(abs(y(3))/sqrt(y(1)^2+y(2)^2));

z=cccgsm(2:4);
zang=atand(abs(z(3))/sqrt(z(1)^2+z(2)^2));

%% Calculate angle
% average b before
bd=bb4n;
da=acosd(sqrt(bd(:,1).^2+bd(:,2).^2)./sqrt(bd(:,1).^2+bd(:,2).^2+bd(:,3).^2))
db=asind(abs(bd(:,3))./sqrt(bd(:,1).^2+bd(:,2).^2+bd(:,3).^2))
dc=atan2(bd(:,3),sqrt(bd(:,1).^2+bd(:,2).^2))*180/pi
%% average angles
da=acosd(sqrt(bb4(:,2).^2+bb4(:,3).^2)./sqrt(bb4(:,2).^2+bb4(:,3).^2+bb4(:,4).^2));
db=asind(abs(bb4(:,4))./sqrt(bb4(:,2).^2+bb4(:,3).^2+bb4(:,4).^2));
dc=atan2(bb4(:,4),sqrt(bb4(:,2).^2+bb4(:,3).^2))*180/pi;
mean(da)
mean(db)
mean(dc)
