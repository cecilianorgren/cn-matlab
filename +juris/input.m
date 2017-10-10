function input(species,n,T,v)
% Input to Juris ES Vlasov code.
%   juris.input(s,n,T,v)
%       Normalised quantities are printed

% Number of species
nSpecies = numel(n);

% Physical constants
e = 1.6022e-19; %units.e;
eps0 = 8.8542e-12; %units.eps0;
mp = 1.6726e-27; %units.mp;
me = 9.1094e-31; %units.me;
c = 299792458; %units.c;

% Normalized density
ntot = sum(n)/2;

% Plasma parameters
for k = 1:nSpecies
    switch species{k}
        case 'e', m(k) = me; q(k) = -1;
        otherwise, m(k) = mp; q(k) = 1;
    end
    op(k) = sqrt(n(k)*1e6*e^2/eps0/m(k)); % rad/s
    vt(k)=sqrt(T(k)*e/m(k)); % root mean square thermal speed
    %vt(k) = c*sqrt(1-1./(T(k).*e./(me*c^2)+1).^2);
    vd(k) = v(k)*vt(1); 
    ld(k) = vt(k)/op(k)/sqrt(2);
    qom(k) = q(k)/(m(k)/me);
end

% Normalization
n_norm = ntot;
v_norm = vt(1);

for k = 1:nSpecies       
    %op(k) = sqrt(n(k)*1e6*e^2/eps0/m(k)); % rad/s
    vt_out(k)=vt(k)/v_norm; % root mean square thermal speed
    %vt(k) = c*sqrt(1-1./(T(k).*e./(me*c^2)+1).^2);
    vd_out(k) = vd(k)/v_norm;
    %ld(k) = vt(k)/op(k)/sqrt(2);
    %qom(k) = q(k)/(m(k)/me);
end





disp(['vdr(): ' num2str(vd_out(1)) ', ' num2str(vd_out(2)) ', ' num2str(vd_out(3))])
disp(['vth(): ' num2str(vt_out(1)) ', ' num2str(vt_out(2)) ', ' num2str(vt_out(3))])
disp(['navg(): ' num2str(n(1)/n_norm) ', ' num2str(n(2)/n_norm) ', ' num2str(n(3)/n_norm)])
disp(['QOM(): ' num2str(qom(1)) ', ' num2str(qom(2)) ', ' num2str(qom(3))])


