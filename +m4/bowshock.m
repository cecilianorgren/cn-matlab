% see Slavin1981
path = '/Users/Cecilia/MATLAB/M4/';
fid = fopen([path 'bowshock_magnetopause_pos_solarwind.txt']);
format ='%*s%*s%f%f%f';
A = textscan(fid,format);
R_BS = sort(A{1}); % RE
R_MP = sort(A{2}); % RE
V_SW = sort(A{3}); % km/s

n = numel(V_SW);

mR_BS = median(R_BS);
mR_MP = median(R_MP);
mV_SW = median(V_SW);
mR_BS_100 = mean(R_BS);
mR_MP_100 = mean(R_MP);
mV_SW_100 = mean(V_SW);
mR_BS_25 = mean(R_BS(1:fix(n/4)));
mR_MP_25 = mean(R_MP(1:fix(n/4)));
mV_SW_25 = mean(V_SW(1:fix(n/4)));
%%
Bz = 0;
Dp=2;
rzero = mR_MP;
alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
        r=rzero*(2./(1+cos(theta))).^alpha;
        xMP=r.*cos(theta);
        yMP=r.*sin(theta);
        ii=find(abs(x)>100);
        xMP(ii)=[];
        yMP(ii)=[];
        
        
%% BS
rstandoff=mR_MP;
        xMP=rstandoff:-0.5:-100;
        rho=sqrt(0.04*(xMP-rstandoff).^2-45.3*(xMP-rstandoff)); % original F/G model adds rstandoff^2=645
        yMP=rho;
        hold(gca,'on')
        plot(gca,xMP,yMP,xBS,yBS,'linewidth',2)
%%
hold(hca,'on')
if 1
    r = @(L,e,theta) L./(1+e*cos(theta));
    theta = (pi/180)*(-90:90);
    e = 1.1;
    L = mR_BS*(1+e);
    xBS = r(L,e,theta).*cos(theta);
    yBS = r(L,e,theta).*sin(theta);
    plot(hca,xBS,yBS)
    %set(gca,'ylim',25*[-1 1],'xlim',[3.5 15])
else
    r = @(mR_BS,e,theta) mR_BS*(1+e)./(1+e*cos(theta));
    theta = (pi/180)*(-90:90);
    e = 1.02;    
    x = r(mR_BS,e,theta).*cos(theta);
    y = r(mR_BS,e,theta).*sin(theta);
    plot(x,y,...
        r(mR_BS_25,e,theta).*cos(theta),r(mR_BS_25,e,theta).*sin(theta)...
        )
    set(gca,'ylim',25*[-1 1],'xlim',[3.5 15])
end
    