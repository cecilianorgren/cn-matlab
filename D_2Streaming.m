% Solving the dispersion relation streaming instabilities

% Frequencies are given in electron backgrund plasma frequency.
% Lengths are given in the inverse electron debye length of the 
% background plasma. Quantities such as vti, vtebeam, opebeam, 
% opi must then be redefined in units of opebg.

% Setting up k (or k*lambda_De_bg)
kmin=0.08; kmax=.4; nk=50;
k=linspace(kmin,kmax,nk);

% Units
mime=1836;

% Plasma components
% Background electrons
te1=1.7;                % keV
ve1=1;                  % vde1/vte1
% Beam electrons
te2=0.37;               % keV
ve2=0.0;                % vde2/vte2
R=10.^(-[ 1 2 3 4]);    % fraction of density in beam
ebs=sqrt(te1*(1-R)./te2./R);  % mult w. k to get the debyelength in in e_bg
% Ions
ti=2.2;                 % eV
vi=0;                   % vdi/vti
is=sqrt(ti*(1-R)*mime/te1); % mult w. k to get the debyelength in in e_bg
K=0;                    % K=1 include ions, K=0, exclude ions

% Functions
zeta=@(x)faddeeva(x)*1i*sqrt(pi);   % Plasma dispersion function
eta=@(x,k,v)( x./k-v );             % Eta
% Dispersion relation
f=@(x,k,v1,v2,v3,ti,te1,te2,R,K,is,ebs)( 1+...
    K./(k.*k.*is.*is)*(1+eta(x,k*is,v2).*zeta(eta(x,k*is,v2)))+... % bg ions
    (1-R)./(k.*k)*(1+eta(x,k,v2).*zeta(eta(x,k,v2)))+... % bg electrons
    R./(k.*k.*ebs.*ebs)*(1+eta(x,k.*ebs,v3).*zeta(eta(x,k.*ebs,v3)))); 
w=[];                               % Omega
omega=zeros(length(R),length(k));

% Solving dispersion relation
for p=1:length(R)
    w=[];                               % Omega
for kk=k
    options=optimset('Display','off'); 
    %x=fsolve(f,1-0.1i,options,kk)*sqrt(2)*kk; 
    x=fsolve(@(x)f(x,kk,vi,ve1,ve2,ti,te1,te2,R(p),K,is(p),ebs(p)),0.01+0.01i,options)*sqrt(2)*kk;
    x=fsolve(@(x)f(x,kk,vi,ve1,ve2,ti,te1,te2,R(p),K,is(p),ebs(p)),x,options)*sqrt(2)*kk;
    
    w=[w,x];
end
    omega(p,:)=w;
end

% Numerical solutions
wre=real(w);
wie=imag(w); 
ore=real(omega);
oie=imag(omega);

% Analytical solution with approximations
%wrt=1.0+1.5.*k.*k; 
%wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3); 
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