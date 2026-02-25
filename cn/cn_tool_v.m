function [x y z corrx corramp] = cn_tool_v(B0,E0,Ni,Ne,gsePos,v,t1,t2,tries)
% Tries different propagation directions 
% and plots phi_B and phi_E.
% Gives back the different vectors x y z, 
% y is propagation direction,
% x is normal direction.

% Define constants
mu0=4*pi*1e-7;
Ne=Ne*1e6;
e=1.6e-19;

% Zoom in
B=cn_toepoch(t1,t2,B0);
E=cn_toepoch(t1,t2,E0);
B=irf_resamp(B,E);

% Filter
[b_av b_hat b_mag]=cn_hat(B);
flh=irf_plasma_calc(b_mag,Ni,0,0,0,'Flh');
flow=0.9*flh;
BAC=irf_filt(B,flow,0,450,5);
EAC=irf_filt(E,flow,0,450,5);

% Set up coordinate systems
theta=1*(0:pi/tries:2*pi-pi/tries)';
size(theta)
x=zeros(tries,3);
y=zeros(tries,3);
z=repmat(b_hat,tries,1);
n_hat=[cos(theta) sin(theta) 0*sin(theta)];
for k=1:tries
    y(k,:)=cn_cross(cn_cross(z(k,:),n_hat(k,:)),z(k,:));
    y(k,:)=y(k,:)/cn_mag(y(k,:));
    x(k,:)=cn_cross(y(k,:),z(k,:));
    x(k,:)=x(k,:)/cn_mag(x(k,:));
end

% Get Phi_B    
scaling=b_mag/e/Ne/mu0/1e18;
Bz=[BAC(:,1) (b_hat*BAC(:,(2:4))')'];
PhiB=[Bz(:,1) Bz(:,2)*scaling];

% Initialize correlations
correlation=zeros(tries,2);


figure;
%h=irf_plot(tries); 
    flag=1;
for k=1:tries   
    % Get E in propagation direction
    Ek=[EAC(:,1) (x(k,:)*EAC(:,(2:4))')'];
    %En=[EAC(:,1) (x(k,:)*EAC(:,(2:4))')'];
    
    % Get Phi_E
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
    PhiE=[PhiE(:,1) PhiE(:,2)-mean(PhiE(:,2))];
       
    % Check correlation
    correlation(k,1)=xcorr(PhiE(:,2),PhiB(:,2),0);
    if flag==1
        sumcorr=sum((PhiE(:,2)-PhiB(:,2)).^2);
        xcorr1=correlation(k,1);
        flag=0;
    end
    correlation(k,2)=sum((PhiE(:,2)-PhiB(:,2)).^2)*abs(xcorr1)/sumcorr;
    
    % Check waveform
    %window=10;
    %[t0x(k),~,~,~,~]=cn_mycorr(PhiE(:,2),PhiB(:,2),window,450);
    %[t0x2(k) corrx2(:,k)]=cn_xcorr(PhiE,PhiB,window,'x');
    
    %irf_plot(h(k),{PhiE,PhiB},'comp')
    %irf_legend(h(k),{'phi_E','phi_B'},[0.90,0.90])
end
i1=find(correlation(:,2)==min(correlation(:,2)))
%il=6
ix=find(correlation(:,1)==max(correlation(:,1)))
Ek=[EAC(:,1) (x(i1,:)*EAC(:,(2:4))')'];
% Get Phi_E
size(PhiE)
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
    PhiE=[PhiE(:,1) PhiE(:,2)-mean(PhiE(:,2))];
h=irf_plot(2);

irf_plot(h(1),{PhiE,PhiB},'comp')
irf_legend(h(1),{'phi_E','phi_B','amp'},[0.90,0.90])
%irf_zoom(h,'x',[PhiE(1,1) PhiE(end,1)])
i2=find(correlation(:,1)==max(correlation(:,1)))
Ek=[EAC(:,1) (x(i2,:)*EAC(:,(2:4))')'];
% Get Phi_E
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
    PhiE=[PhiE(:,1) PhiE(:,2)-mean(PhiE(:,2))];
    irf_plot(h(2),{PhiE,PhiB},'comp')
irf_legend(h(2),{'phi_E','phi_B','xcorr'},[0.90,0.90])

correlation;
corrx=correlation(:,1);
corramp=correlation(:,2);