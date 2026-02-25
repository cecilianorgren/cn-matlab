function input(varargin)
% Prepare input to xes1 runs.
% Time normalizations are with respect to total electron density and
% velocity normalizations with respect to first species.
%   xes.input(n,T,vd,species)
%       n - cc
%       T - eV
%       vd - km/s       
%       species - 'e' or 'i'/'p'
%
Units=irf_units; % read in standard units
normspecies = 1;

Me=Units.me;
Mp=Units.mp;
c=Units.c;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;

n = varargin{1}*1e6; % m^-3
T = varargin{2}; % eV
vd = varargin{3}*1000; % m/s
species = varargin{4};

nSpecies = numel(species);
ntot = sum(n)/2;
wpnorm =  sqrt(ntot*e^2/Me/epso); % rad/s    

for k = 1:nSpecies
    if any([strcmp(species{k},'i') strcmp(species{k},'p')]); % ions
        m(k) = Mp;
    else
        m(k) = Me;
    end
    wp(k) =  sqrt(n(k)*e^2/m(k)/epso); % rad/s    
    vt(k) = c*sqrt(1-1./(T(k).*e./(m(k)*c^2)+1).^2)/sqrt(2);              % m/s (relativ. correct), particle with Vte has energy e*Te        
end        
for k = 1:nSpecies
    disp(['wp' num2str(k) ' = ' num2str(wp(k)/wpnorm,'%1.3f') ', vt' num2str(k) ' = (k_B*T_s/m_s)^1/2 = ' num2str(vt(k)/vt(normspecies),'%1.3f') ', vd' num2str(k) ' = ' num2str(vd(k)/vt(normspecies),'%1.3f')])    
end

%disp(vt/1000)
%disp(vd/1000)