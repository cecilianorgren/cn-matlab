function Ln = cn_ln(t1,t2,ExB,flh,PotNe,peaNe,M)
% c_eval('ExB?=cn_toepoch(t1,t2,ExB?);',3:4);
dt=ExB(2,1)-ExB(1,1);
fs=1/dt;
flow=flh*0.4;
fhigh=180;
% Taking a slightly larger interval to make the filtering
t11=t1;
t11(6)=t11(6)-1;
t22=t2;
t22(6)=t22(6)+1;

ExB=cn_toepoch(t11,t22,ExB);
ExBDC=irf_filt(ExB,0,flow,fs,3);
ExBDC=cn_toepoch(t1,t2,ExBDC);

bExB=M*ExBDC(:,2:4)';
bExB=[ExBDC(:,1) bExB'];

% Filter away all small density fluctuations
%PotNe=irf_filt(PotNe,0,flow,fs,3);
%PotNe=[PotNe(:,1) PotNe(:,2)/20];
%figure;irf_plot(PNe)

%figure;irf_plot(bExB)
dx=dt*sum(bExB(:,2),1); % y-component (normal) since time is taken away
n1=cn_toepoch(t1,PotNe);n1(2);
n2=cn_toepoch(t2,PotNe);n2(2);
dn=n2-n1;
dn=dn(2);
%n=cn_toepoch(t1,peaNe);
%n=n(2);
n_all=cn_toepoch(t1,t2,PotNe);
n_all(:,2);
n=sum(n_all(:,2))/size(n_all,1); % 20 is from calibrating with PEACE
size(n_all,1);
size(n_all);


Ln=abs(n*dx/dn);
