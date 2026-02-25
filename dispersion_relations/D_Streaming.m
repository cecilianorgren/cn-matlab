
% Setting up k (or k*lambda_De)
kmin=0.05;                             
kmax=0.4;                         
k=linspace(kmin,kmax,100);

% Units
mime=1836;
% Things to define
te1=1.7;            % keV
te2=0.17;           % keV
ti=2;
vi=0;               % vdi/vti
ve1=0;              % vde1/vte1
ve2=2.6;           % vde2/vte2
R=0.1;
K=0;                % K=1 include ions, K=0, exclude ions

zeta=@(x)faddeeva(x)*1i*sqrt(pi);   % Plasma dispersion function
eta=@(x,k,v)( x./k-v );             % Eta
% Dispersion relation
f=@(x,k,v1,v2,v3,ti,te1,te2,R)( 1+K./(k*k)*(1+eta(x,k*sqrt(mime*te1/ti),v1)*zeta(eta(x,k*sqrt(mime*te1/ti),v1)))+...
                   (1-R)./(k*k)*(1+eta(x,k,v2)*zeta(eta(x,k,v2)))+...
                       R./(k*k)*(1+eta(x,k/sqrt((1-R)*te2/R/te1),v3)*zeta(eta(x,k*sqrt((1-R)*te2/R/te1),v3)))); 
w=[];                               % Omega


% Solving dispersion relation
for kk=k
    options=optimset('Display','off'); 
    %x=fsolve(f,1-0.1i,options,kk)*sqrt(2)*kk; 
    x=fsolve(@(x)f(x,kk,vi,ve1,ve2,ti,te1,te2,R),0.1-0.01i,options)*sqrt(2)*kk;
    w=[w,x];
end

% Numerical solutions
wre=real(w);
wie=imag(w); 

% Analytical solution with approximations
%wrt=1.0+1.5.*k.*k; 
%wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3); 
plot(k,wre,'-b',k,-wie,'-g'); 
%legend('\omega_r','-\gamma','\omega_r','-\gamma','Location','SouthEas t');
grid on; 
title(strcat('Three stream instability '...
    ,10,'numerical computation'));
xlabel('k\lambda_D');
ylabel('\omega/\omega_p'); 
%xlim([kmin,kmax]);  
%ylim([0.0001,100]);