clear all
close all

% simulation box length
L=2*2*pi;
% time step
DT=0.1;
% number of computational cycles
NT=150;
% number of grid points
NG=64;
% grid spacing
dx=L/NG;
% number of particles
N=300000;
% plasma frequency
WP=1;
% electron charge to mass ratio
QM=-1;
% beam velocity
V0=0.0;
% thermal velocity
VT=1.0;
% perturbation amplitude and mode
XP1=0.1; 
V1=0.0;
mode=1;
% computational particle charge
Q=WP^2/(QM*N/L);
%background charge given by background (not moving) ions
rho_back=-Q*N/L;

% 2 Stream instability set-up
% uniform particle distribution in space
xp=linspace(0,L-L/N,N)';
% Maxwellian velocity 
vp=VT*randn(N,1);
pm=[1:N]';
pm=1-2*mod(pm,2);
% add the beam drift velocity
vp=vp+pm.*V0;

% Perturbation in particle positions
xp=xp+0.1*sin(2*pi*xp/L*mode);

% auxiliary vectors and matrices
p=1:N;p=[p p];
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);

% history vectors aand matrices
histTotEnergy = [];
histKinEnergy = [];
histPotEnergy = [];
histMomentum = [];
histThermalVelocity = [];
histSpectrum = [];
histEfield = [];
histPhi    = [];
histPARTx = [];
histPARTv = [];



for it=1:NT
   
   % update particle position xp
   xp=xp+vp*DT;  
   % apply the boundary conditions
   out=(xp<0); xp(out)=xp(out)+L;
   out=(xp>=L);xp(out)=xp(out)-L;

   % interpolation from particle to grid
   g1=floor(xp/dx-.5)+1;
   g=[g1;g1+1];
   fraz1=1-abs(xp/dx-g1+.5);
   fraz=[fraz1;1-fraz1]	;
   out=(g<1);g(out)=g(out)+NG;
   out=(g>NG);g(out)=g(out)-NG;
   mat=sparse(p,g,fraz,N,NG);
   % calculate the charge density
   rho=full((Q/dx)*sum(mat))'+rho_back;

   % calculate the electrostatic potential'
   Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
   
   % calculate the electric field from the electrostatic potential
   Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
   
   % interpolation grid -> particle and velocity update
   vp=vp+mat*QM*Eg*DT;
   
   % Calculate different kind of energies
   Ekin   = 0.5*abs(Q)*sum(vp.^2);
   % Potential energy
   Efield = 0.5*sum(Eg.^2)*dx;
   % Total energy
   Etot   =  Ekin + Efield ;
   
   % Total momentum
   histMomentum  = [histMomentum  sum(abs(Q)*vp)];
   
   % history of various energies
   histKinEnergy = [histKinEnergy Ekin];
   histPotEnergy   = [histPotEnergy Efield];
   histTotEnergy = [histTotEnergy Etot];
   % Take only the even particle (just one beam)
   histThermalVelocity = [histThermalVelocity cov(vp(1:2:end))];
   histEfield = [histEfield Eg];
   histPhi    = [histPhi Phi];
   
   % take the Fourier transform of the electric field on the grid
   NFFT = 2^nextpow2(length(Eg)); % Next power of 2 from length of Eg
   Y = fft(Eg,NFFT)/length(Eg);
   histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];
   
   % take the phase space after a given iteration (taking every 10th particle)
   if (it > 500)
             histPARTx = [histPARTx; xp(1:10:end)'];
			 histPARTv = [histPARTv; vp(1:10:end)'];
   endif
   
end

% useful varaibles for plotting
time = linspace(0,NT*DT,NT);
k = 2*pi*(1/(2*dx))*linspace(0,1,NFFT/2+1);
space = linspace(0,L,NG);
