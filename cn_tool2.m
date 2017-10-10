function [x y z] = cn_tool2(B,E,Ni,Ne,gsePos,v,t1,t2,tries)
% Define constants
mu0=4*pi*1e-7;
Ne=Ne*1e6;
e=1.6e-19;

% Zoom in
B=cn_toepoch(t1,t2,B);
E=cn_toepoch(t1,t2,E);

% Filter
[b_av b_hat b_mag]=cn_hat(B);
flh=irf_plasma_calc(b_mag,Ni,0,0,0,'Flh');
flow=0.5*flh;
BAC=irf_filt(B,flow,0,450,5);
EAC=irf_filt(E,flow,0,450,5);

% Set up coordinate systems
theta=(0:pi/tries:pi-pi/tries)';
x=zeros(tries,3);
y=zeros(tries,3);
z=repmat(b_hat,tries,1);
n_hat=[cos(theta) sin(theta) 0*sin(theta)];
for k=1:tries
    x(k,:)=cn_cross(cn_cross(z(k,:),n_hat(k,:)),z(k,:));
    x(k,:)=x(k,:)/cn_mag(x(k,:));
    y(k,:)=cn_cross(z(k,:),x(k,:));
    y(k,:)=y(k,:)/cn_mag(y(k,:));
end

% Get Phi_B    
scaling=b_mag/e/Ne/mu0/1e18;
Bz=[BAC(:,1) (b_hat*BAC(:,(2:4))')'];
PhiB=[Bz(:,1) Bz(:,2)*scaling];

h=irf_plot(tries);
for k=1:tries    
    % Get E in propagation direction
    Ek=[EAC(:,1) (y(k,:)*EAC(:,(2:4))')'];
    En=[EAC(:,1) (x(k,:)*EAC(:,(2:4))')'];
    
    % Get Phi_E
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -v*PhiE(:,2)];
    
    % Get dt
    %Pos=cn_toepoch(t1,gsePos);
    %Posy=y*Pos';
    
    % Check waveform
    %window=10;
    %[t0x(k),~,~,~,~]=cn_mycorr(PhiE(:,2),PhiB(:,2),window,450);
    %[t0x2(k) corrx2(:,k)]=cn_xcorr(PhiE,PhiB,window,'x');
    
    irf_plot(h(k),{PhiE,PhiB},'comp')
    irf_legend(h(k),{'phi_E','phi_B'},[0.90,0.90])
end

