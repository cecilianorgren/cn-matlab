% input data
Te1 = 1600;
Te2 = [60];
Ti = [1600];
S = 0.3:0.1:1.4;
R = [0.4];
k = 0.01:0.02:0.5; nk = numel(k);
n = 0.06;

e = 1.6022e-19; %units.e;
eps0 = 8.8542e-12; %units.eps0;
mp = 1.6726e-27; %units.mp;
me = 9.1094e-31; %units.me;
c = 299792458; %units.c;

omega_pi = sqrt(n*1e6*e^2/eps0/mp); % rad/s
omega_pe = sqrt(n*1e6*e^2/eps0/me); % rad/s
omega_pe1 = sqrt(n*1e6*e^2/eps0/me)*sqrt(1-R); % rad/s
omega_pe2 = sqrt(n*1e6*e^2/eps0/me)*sqrt(R); % rad/s

vthe1 = c*sqrt(1-1./(Te1.*e./(me*c^2)+1).^2);
vthe2 = c*sqrt(1-1./(Te2.*e./(me*c^2)+1).^2);
vthi = c*sqrt(1-1./(Ti.*e./(mp*c^2)+1).^2);

vde1 = 0; % m/s
vde2 = S*vthe1; % m/s
vdi = 0; % m/s

lamD1 = (vthe1')*(1./omega_pe1)/sqrt(2); % m
lamD2 = (vthe2')*(1./omega_pe2)/sqrt(2); % m

lamD = (vthe1')*(1./omega_pe)/sqrt(2); % m
%l_norm = v_norm/w_norm/sqrt(2);
% omega_bune = omega_pe^(1/3)*omega_pi^(2/3)*16^(-1/3);

normal_frequency = omega_pe; % total electron plasma frequency
normal_length = lamD; % Vte./Wpe/sqrt(2); % total density, background temperature
normal_frequency_string = '\omega_{pe}';
normal_length_string = '\lambda_D';
% run the solver

TolFun = 1e-8;
%TolFun = 1e-12;
tic
nw = 100;
iR=1;
iTi=1;
iTe1=1;
iTe2=1;
iv = 1;
fv = [omega_pi,omega_pe*sqrt(1-R(iR)),omega_pe*sqrt(R(iR)),vthi(iTi),vthe1,vthe2(iTe2),vdi,vde1,vthe1*S(iv),lamD];

wr0 = [0];
wi0 = [0];
wi = nan(1,nk);
wr = nan(1,nk);
for ik = 1:nk
    wrvec = wr0(ik)+linspace(-200,200,nw);
    wivec = wi0(ik)+linspace(-200,200,nw);
    [WI,WR] = meshgrid(wivec,wrvec); 
    WR(WR<0) = NaN;
    for pp = 1:nw;
        for nn = 1:nw;           
            w_res = csolve_disprel([WI(pp,nn) WR(pp,nn)],k(ik),fv);
            WRres(pp,nn) = real(w_res);
            WIres(pp,nn) = imag(w_res);
        end
    end    
    ressum = abs(WRres)+abs(WIres);
    minreswr = WR(find(ressum==min(ressum(:))));
    minreswi = WI(find(ressum==min(ressum(:))));
    disp(['residual = [' num2str(minreswr) ' + ' num2str(minreswi) 'i]'])
    if isempty(minreswr)
        wr(ik)=wr(ik-1);
        wr(ik)=wr(ik-1);
    else
        wr(ik)=WR(find(ressum==min(ressum(:))));
        wi(ik)=WI(find(ressum==min(ressum(:))));
    end
    if ik == 1 % interpolate from zero.
        dwr = wr(ik)-0;
        dwi = wi(ik)-0;
    else   
        dwr = wr(ik)-wr(ik-1);
        dwi = wi(ik)-wi(ik-1);
    end    
    wr0(ik+1) = wr(ik)+dwr;
    wi0(ik+1) = wi(ik)+dwi;
end


subplot(1,2,1)
plot(k,wr/normal_frequency,k,wi/normal_frequency)
legend('w_r','w_i')
subplot(1,2,2)
plot(k,wr,k,wi)
legend('w_r','w_i')
%% save
saveName = [datestr(now,'yyyy-mm-ddTHHMMSS') '_beam_solver_' num2str(numel(k)) 'kx' num2str(numel(S)) 'Sx' num2str(numel(R)) 'Rx' num2str(numel(Te2)) 'Te2x' num2str(numel(Ti)) 'Ti'];
save(['/home/cecilia/Research/BeamSolver/' saveName])
