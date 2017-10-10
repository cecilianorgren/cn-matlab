% Solving the dispersion relation streaming instabilities

% Frequencies are given in electron backgrund plasma frequency.
% Lengths are given in the inverse electron debye length of the 
% background plasma. Quantities such as vti, vtebeam, opebeam, 
% opi must then be redefined in units of opebg.

% Clearing old variables
clear R;

% Setting up k (or k*lambda_De_bg)
kmin=0.08; kmax=.8; nk=20;
k=linspace(kmin,kmax,nk);

% Units
mime=1836;
memi=1/1836;

% Plasma components
% 1. Background electrons
te1=1.7;                % keV
ve1=0.1;                  % vde1/vte1

% 2. Beam electrons
te2=0.37;               % keV
ve2=6.92;                % vde2/vte2
R=10.^(-[0 0 0 0]);    % fraction of density in beam

% Ions
ti=2.2;                 % eV
vi=0.1;                   % vdi/vti
K=1;                    % K=1 include ions, K=0, exclude ions

% Rescaling factors


% Functions
zeta=@(x)(faddeeva(x)*1i*sqrt(pi));   % Plasma dispersion function
eta=@(x,k,v,s)( x.*s./k-v );             % Eta
t=@(R,t1,t2)(R*t2+(1-R)*t1);                  % Average temperature
si=sqrt(memi*(t(R,te1,te2)/ti));
se1=sqrt(t(R,te1,te2)/te1);
se2=sqrt(t(R,te1,te2)/te2);

chi_i=@(o,k,K,v,s)( (2*K*1836*s^2/(k.^2).*(1+eta(o,k,v,s)*zeta(eta(o,k,v,s)))) );
chi_e1=@(o,k,R,v,s)( (2*(1-R)*s^2/(k.^2).*(1+eta(o,k,v,s)*zeta(eta(o,k,v,s)))) );
chi_e2=@(o,k,R,v,s)( (2*R*s^2/(k.^2).*(1+eta(o,k,v,s)*zeta(eta(o,k,v,s)))) );
% Dispersion relation
%f=@(x,k,v1,v2,v3,ti,te1,te2,R,K,is,ebs)( 1+...
%    K.*is.*is./(k.*k)*(1+eta(x,k*is,v2).*zeta(eta(x,k*is,v2)))+... % bg ions
%    (1-R)./(k.*k)*(1+eta(x,k,v2).*zeta(eta(x,k,v2)))+... % bg electrons
%    R./(k.*k.*ebs.*ebs)*(1+eta(x,k.*ebs,v3).*zeta(eta(x,k.*ebs,v3)))); 
D=@(o,k,R,K,vi,ve1,ve2,si,se1,se2)...
    (1+chi_i(o,k,K,vi,si)+chi_e1(o,k,R,ve1,se1)+chi_e2(o,k,R,ve2,se2));
w=[];                               % Omega
omega=zeros(length(R),length(k));

% Solving dispersion relation
for p=1:length(R)
for q=1:nk
    options=optimset('Display','off'); 
    x=fsolve(@(x)D(x,k(q),R,K,vi,ve1,ve2,si(p),se1(p),se2(p)),0.1-0.1i,options);
    omega(p,q)=fsolve(@(x)D(x,k(q),R,K,vi,ve1,ve2,si(p),se1(p),se2(p)),0.01-0.01i,options);
end
end

% Numerical solutions
ore=real(omega);
oie=imag(omega);

% Analytical solution with approximations
%wrt=1.0+1.5.*k.*k; 
%wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3); 
figure(1);
h1=subplot(2,1,1);
plot(h1,k,oie(1,:),'-b',k,oie(2,:),'-g',k,oie(3,:),'-r',k,oie(4,:),'-k'); 
%legend('\omega_r','-\gamma','\omega_r','-\gamma','Location','SouthEas t');
grid on; 
title(strcat('Three stream instability '...
    ,10,'numerical computation'));
xlabel(h1,'k\lambda_D');
ylabel(h1,'\gamma/\omega_p'); 
xlim([kmin,kmax]);  
%ylim([-0.1,0.1]);

h2=subplot(2,1,2);
plot(h2,k,ore(1,:),'-b',k,ore(2,:),'-g',k,ore(3,:),'-r',k,ore(4,:),'-k'); 
xlabel(h2,'k\lambda_D');
ylabel(h2,'\omega/\omega_p'); 
xlim([kmin,kmax]);  
%ylim([-0.1,0.1]);

linkaxes([h1 h2],'x')