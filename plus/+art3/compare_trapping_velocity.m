cd /Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/
doLoad = 1;
doClearup = 1;

if doLoad
    %load('Ti-2000eV_eventhickergrid')
    %load('Ti-500eV_eventhickergrid')
    load('Te1-2000_Te2-120_Ti-2000eV_R15S15k40')
    nR = numel(R);
    nv = numel(S);
    nk = numel(k)

    B = 25;
    n = 0.06;
    no = 0;

    % Ions
    omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
    vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
    vdi = 0;
    % Total electron plasma frequency
    omega_pe = irf_plasma_calc(B,n,no,Te1,Ti,'Fpe')*2*pi; % rad/s
    % n1= n*(1-R); % background
    % n2= n*R;     % beam 
    % Background electrons
    % omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
    vthe1 = irf_plasma_calc(B,n,no,Te1,Ti,'Vte'); % m/s 
    vde1 = 0; % m/s
    % Beam electrons
    % omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
    vthe2 = irf_plasma_calc(B,n,no,Te2,Ti,'Vte'); % m/s 
    vde2 = S*vthe1;
    % Other
    lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m
    omega_bune = omega_pe^(1/3)*omega_pi^(2/3)*16^(-1/3);

    normal_frequency = omega_pe; % total electron plasma frequency
    normal_length = lamD; % Vte./Wpe/sqrt(2); % total density, background temperature
    normal_frequency_string = '\omega_{pe}';
    normal_length_string = '\lambda_D';
end
if doClearup
    %% clear up the matrix from phoney values
    k_store = repmat(k',1,nv,nR);
    vph_store = x_real_store/normal_frequency./k_store/sqrt(2);
    x_imag_store(x_real_store<0) = NaN;
    vph_store(x_real_store<0) = NaN;
    x_real_store(x_real_store<0) = NaN;
    x_real_store(x_real_store/normal_frequency>2.4) = NaN;
    %% find the point where the growth is max
    % find the max for each R and S, i.e. one per k-vector
    x_imag_max_store_3D = zeros(size(x_imag_store));
    x_imag_max_store_2D = zeros(nv,nR);
    x_imag_max_store_2Dvalues = zeros(nv,nR);
    vph_max_store = zeros(nv,nR);
    k_max_store  = zeros(nv,nR);
    %(x_imag_max_store_2D(iv,:)
    for iR = 1:nR
        for iS = 1:nv
            [ind,~] = find(x_imag_store(:,iS,iR)==max(x_imag_store(:,iS,iR)));

            x_imag_max_store_3D(ind,iS,iR)=1;    
            if ~isempty(ind)
                x_imag_max_store_2D(iS,iR)=ind;
                vph_max_store(iS,iR) = vph_store(ind,iS,iR);
                k_max_store(iS,iR) = k_store(ind,iS,iR);
                x_imag_max_store_2Dvalues(iS,iR) = x_imag_store(ind,iS,iR);
            else
                x_imag_max_store_2D(iS,iR)=NaN;
                vph_max_store(iS,iR) = NaN;
                k_max_store(iS,iR) = NaN;
                x_imag_max_store_2Dvalues(iS,iR) = NaN;
            end
            disp(['iR=' num2str(iR) ', iS=' num2str(iS) ', max_ind=' num2str(ind) ', vph_max_ind=' num2str(vph_max_store(iS,iR))])        
        end
    end
end   
vph_plot = vph_max_store;
vph_plot(vph_max_store>2) = NaN;
vph_plot(vph_max_store<1e-3) = NaN;


f_tot = 1;
f_vthebeam = 0;
units = irf_units;
vthebeam = cn_eV2v(Te2,'eV'); % km/s
vthebg = cn_eV2v(Te1,'eV'); % km/s
vbeam = repmat(S',1,numel(R))*vthebg;

vphi = f_tot*(vbeam-f_vthebeam*vthebeam-vph_plot*vthebg);
% Trapping velocity from Omura1996
Rmat = repmat(R,numel(S),1);

a1 = vbeam-1*vthebeam;
%a1 = vbeam-vph_plot*vthebg;
a2 = 8/3*Rmat;
a3 = 9*vbeam*vthebeam/2./Rmat./a1./a1;
vT = a1*sqrt(a2*(sqrt(1+a3)-1))*1e-3;
%vT = (vbeam-1*vthebeam)*(sqrt(8/3*Rmat*sqrt(1+9*vbeam*vthebeam/2./Rmat/(vbeam-vthebeam)/(vbeam-vthebeam))-1))
% vphi = vT;
vph = vph_plot*vthebg;


nrows=2;
ncols=2;
subplot(nrows,ncols,1)
plot(R,vT,'o')
xlabel('R')
ylabel('v_T')

subplot(nrows,ncols,2)
plot(S,vT','o')
xlabel('S')
ylabel('v_T')

subplot(nrows,ncols,3)
plot(R,vphi,'o')
xlabel('R')
ylabel('v_{\phi}')

subplot(nrows,ncols,4)
plot(S,vphi','o')
xlabel('S')
ylabel('v_{\phi}')

