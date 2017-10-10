% D_3cold.m
% Solves cold plasma "3"-stream instability
clear wi wr w k R
% Setting up k (or k*lambda_De_bg)
kmin=0.01; kmax=0.8; nk=50;
k=linspace(kmin,kmax,nk);

% Units
mime=1836;

% Plasma components
% Background electrons
tebg=1.7;%1.7;           % keV
vebg=1;                  % vde1/vte1
% Beam electrons
teb=0.37;               % keV
veb=8;                % vde2/vte2
R=10.^(-[0 1 2 3 ]);    % fraction of density in beam
% Ions
ti=2.2;                 % eV
vi=1;                   % vdi/vti
K=1;                    % K=1 include ions, K=0, exclude ions

% "Average" T
te=@(tebg,teb,R)(R*teb+(1-R)*tebg);
si=@(v,ti,teb,tebg,R)(v*sqrt(ti./te(tebg,teb,R)/1836));
seb=@(v,teb,tebg,R)(v*sqrt(teb./te(tebg,teb,R)));
sebg=@(v,teb,tebg,R)(v*sqrt(tebg./te(tebg,teb,R)));

% Dispersion relation
D=@(x,k,vi,vebg,veb,ti,teb,tebg,R,K)(1-...
    K*(1/1836)./(x-k.*si(vi,ti,teb,tebg,R)).^2-...
    R/(x-k*seb(veb,teb,tebg,R)).^2-...
    (1-R)/(x-k*sebg(vebg,teb,tebg,R)).^2);

w=zeros(nk,length(R));                % Omega

% Solving dispersion relation
for p=1:length(R) 
for q=1:nk
    options=optimset('Display','off'); 
    %x=fsolve(f,1-0.1i,options,kk)*sqrt(2)*kk; 
    x=fsolve(@(x)D(x,k(q),vi,veb,vebg,ti,teb,tebg,R(p),K),.1+.1i,options);
    w(q,p)=fsolve(@(x)D(x,k(q),vi,veb,vebg,ti,teb,tebg,R(p),K),x,options);
end
end

wi=imag(w);
wr=real(w);

figure(1);
color={'b','g','r','k','c','y'};
h1=subplot(2,1,1);
h2=subplot(2,1,2);
for p=1:length(R)
    plot(h1,k,wi(:,p),char(color{p})); hold(h1,'on'); 
    plot(h2,k,wr(:,p),char(color{p})); hold(h2,'on'); 
end
hold(h1,'off');hold(h2,'off')
%legend('\omega_r','-\gamma','\omega_r','-\gamma','Location','SouthEas t');
grid on; 
title(strcat('Three stream instability '...
    ,10,'numerical computation'));
xlabel(h1,'k\lambda_D');
ylabel(h1,'\gamma/\omega_{pe}'); 
xlabel(h2,'k\lambda_D');
ylabel(h2,'\omega/\omega_{pe}'); 
xlim([kmin,kmax]); 
