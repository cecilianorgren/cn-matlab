function [k v dt n] = cn_tool(B,E,Ni,Ne,t1,t2,tries)
mu0=4*pi*1e-7;
Ne=Ne*1e6;
e=1.6e-19;

k=1; v=1; dt=1;
[b_av b_hat b_mag]=cn_hat(B);

flh=irf_plasma_calc(b_mag,Ni,0,0,0,'Flh')
flow=0.5*flh;

B=cn_toepoch(t1,t2,B);
E=cn_toepoch(t1,t2,E);
BAC=irf_filt(B,flow,0,450,5);
EAC=irf_filt(E,flow,0,450,5);


z=b_hat;
theta=0:2*pi/tries:(2*pi-2*pi/tries);
theta=theta';

% Get Phi_B    
Bz=[BAC(:,1) (z*BAC(:,(2:4))')'];
PhiB=[Bz(:,1) Bz(:,2)*b_mag/e/Ne/mu0/1e18];

h=irf_plot(tries)
isub=1;

for k=1:tries
    
    
end


n=[cos(theta) sin(theta) 0*sin(theta)];
for k=1:tries
    n_hat=[cos(theta(k)) sin(theta(k)) 0*sin(theta(k))];
    x=cn_cross(cn_cross(z,n_hat),z);
    x=x/cn_mag(x)
    y=cn_cross(z,x);
    y=y/cn_mag(y)
    M=[x;y;z];
    
    % Get E in propagation direction
    Ek=[EAC(:,1) (y*EAC(:,(2:4))')'];
    En=[EAC(:,1) (x*EAC(:,(2:4))')'];
    %irf_plot({Ek,Ek},'comp')
    
    % Get Phi_E
    PhiE=irf_integrate(Ek);
    PhiE=[PhiE(:,1) -PhiE(:,2)];
    
    
    
    % Check waveform
    window=10;
    [t0x(k),~,~,~,~]=cn_mycorr(PhiE(:,2),PhiB(:,2),window,450);
    [t0x2(k) corrx2(:,k)]=cn_xcorr(PhiE,PhiB,window,'x');
    
    irf_plot(h(k),{PhiE,Bz},'comp')
end


for k=1:tries
    
end

n_hat=cn_m_trans([1 n_hat],[x;y;z],-1);
dt=t0x';
k=t0x2';
