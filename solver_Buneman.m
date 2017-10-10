% Two-stream script

% Range of parameters
nk = 60;
nv = 10;
k = linspace(0.3,1,nk);
S = linspace(.4,1.1,nv);
S = [1.37];
nv = numel(S);
R = .3;
nR = numel(R);
incIons=1;

% Physical variables:
B = 25;
n = 0.06;

n1= n*(1-R);
n2= n*R;
no = 0;
Te1 = 1600; TeK1 = Te1*11604; % eV -> K
Te2 = 1;   TeK2 = Te2*11604; % eV -> K
Ti = 2200;  TiK  = Ti*11604; % eV -> K

% Ions
omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
vdi = 0;
% Total electrons
omega_pe = irf_plasma_calc(B,n1+n2,no,Te1,Ti,'Fpe')*2*pi; % rad/s
% Background electrons
omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
vthe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Vte'); % m/s 
vde1 = 0; % m/s
% Beam electrons
omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
vthe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Vte'); % m/s 
vde2 = S*vthe1;
% Other
lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m
omega_bune = omega_pe1^(1/3)*omega_pi^(2/3)*16^(-1/3);

normal_frequency = omega_pe;
normal_length = lamD;
normal_frequency_string = '\omega_{pe}';
normal_length_string = '\lambda_D';

% set up figure for drawing results
doPlot=1;
if doPlot 
    figure(10)    
    for nPlot = 1:2     
        h(nPlot) = subplot(1,2,nPlot); hold(h(nPlot),'on');        
        xlabel(h(nPlot),['k' normal_length_string],'Fontsize',14);
        ylabel(h(nPlot),['\omega/'  normal_frequency_string],'Fontsize',14);    
        set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
    end
    title(h(1),'\Re (\omega)')
    title(h(2),'\Im (\omega)')
end


nTryLow = 1;
nTryHigh = 0;

rights = zeros(1,nk);

% Starting guesses, can find it from solver_BunemanRandom
x_real_store = ones(1,nv)*0.01*normal_frequency;
x_imag_store = ones(1,nv)*0.05*normal_frequency;

tic;
for ik=1:nk; %disp(['k=' num2str(k(ik))])        
    for iv=1:nv % For each k*lD
        for iR = 1:nR % Different beam densities
            fv = [incIons*omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2(iv),lamD];
            af = @(temp) myfun_20070831_EH(temp,k(ik),fv);%,lamD,omega_pi,vthi,omega_pe,);
            [x,FVAL,EXITFLAG] = fsolve(af,[x_real_store(end,iv) x_imag_store(end,iv)],optimoptions('fsolve','MaxIter',1000,'TolFun',1e-13,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));
            %if(x(1)>0 && x(1)<normal_frequency*0.8)
                x_real_store(ik,iv) = x(1);
                x_imag_store(ik,iv) = x(2);
            %end
            if doPlot
                if(EXITFLAG==1) % Found something, now plot            
                    plot(h(1),k(ik),x(1)/normal_frequency,'g+','Markersize',3);             
                    plot(h(2),k(ik),x(2)/normal_frequency,'g+','Markersize',3); hold on;
                    % only plot for reasonably omega_i > -0.6            

                    if(x(2)/normal_frequency>-0.6)
                        if ( abs(x(1)/normal_frequency)<0.8)                    
                            plot(h(1),k(ik),x(1)/normal_frequency,'r+','Markersize',6); hold on; 
                            plot(h(2),k(ik),x(2)/normal_frequency,'r+','Markersize',6); hold on;                
                        elseif 0;( abs(x(1)/normal_frequency)>0.8) % before it was 0.1                    
                            plot(h(1),k(ik),x(1)/normal_frequency,'bd','Markersize',6); hold on;                    
                            plot(h(2),k(ik),x(2)/normal_frequency,'bd','Markersize',6); hold on;    
                        end
                        if (x(2)/normal_frequency>0.0001)                    
                            rights(ik)=rights(ik)+1;
                            plot(h(1),k(ik),x(1)/normal_frequency,'ch','Markersize',6); hold on;                    
                            plot(h(2),k(ik),x(2)/normal_frequency,'ch','Markersize',6); hold on; 
                            disp(['x_real=' num2str(x_real_store(end,iv),'%.4f') ', x_imag=' num2str(x_imag_store(end,iv),'%.4f')])
                        end
                    end                   
                    drawnow  
                end                
            end
        end % end nR
    end % end nv
    rights;
end % end nk
toc;

if 0
    set(h(1),'ylim',[0 1],'xlim',[0 2]);
    set(h(2),'ylim',[0 0.5],'xlim',[0 2])
end
