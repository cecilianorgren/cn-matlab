% Solving the dispersion relation streaming instabilities

% Frequencies are given in electron backgrund plasma frequency.
% Lengths are given in the inverse electron debye length of the 
% background plasma. Quantities such as vti, vtebeam, opebeam, 
% opi must then be redefined in units of opebg.

% Clearing variables
clear R
% Setting up k (or k*lambda_De_bg)
kmin=0.08; kmax=.4; nk=30;
k=linspace(kmin,kmax,nk);

% Units
mime=1836;

% Plasma componentsThings to define
% Background electrons
te1=1.7;%1.7;           % keV
ve1=.5;                  % vde1/vte1
% Beam electrons
te2=0.37;               % keV
ve2=-6.92;                % vde2/vte2
R=10.^(-[ 1 2 3 4]);    % fraction of density in beam
% Ions
ti=2.2;                 % eV
vi=0;                   % vdi/vti
K=0;                    % K=1 include ions, K=0, exclude ions



zeta=@(x)faddeeva(x)*1i*sqrt(pi);   % Plasma dispersion function
eta=@(x,k,v,s)( x.*s./k-v );             % Eta
% Dispersion relation
f=@(x,k,v1,v2,v3,ti,te1,te2,R,K)( 1+...
    2*K.*(te1/ti./(k.*k))*(1+eta(x,k,v2,sqrt(te1/ti/1836)).*zeta(eta(x,k,v2,sqrt(te1/ti/1836))))+... % bg ions
    2*(1-R)./(k.*k)*(1+eta(x,k,v2,1).*zeta(eta(x,k,v2,1)))+... % bg electrons
    2*R.*(te1/te2)./(k.*k)*(1+eta(x,k,v3,sqrt(te1/te2)).*zeta(eta(x,k,v3,sqrt(te1/te2))))); 

omega=zeros(length(R),length(k));

% Solving dispersion relation
for p=1:length(R)
for q=1:nk
    options=optimset('Display','off'); 
    %x=fsolve(f,1-0.1i,options,kk)*sqrt(2)*kk; 
    x=fsolve(@(x)f(x,k(q),vi,ve1,ve2,ti,te1,te2,R(p),K),0.01-0.01i,options)*sqrt(2)*kk;
    omega(p,q)=fsolve(@(x)f(x,k(q),vi,ve1,ve2,ti,te1,te2,R(p),K),x,options)*sqrt(2)*kk;
end
end

% Numerical solutions
wre=real(w);
wie=imag(w); 
ore=real(omega);
oie=imag(omega);

% Analytical solution with approximations
%wrt=1.0+1.5.*k.*k; 
%wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3); 
figure(2);
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