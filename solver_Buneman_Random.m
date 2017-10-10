% Two-stream script
% Variable parameters

% Physical variables:
vde1=0;
event = 111;
no = 0;
switch event
    case 1
        S=1.6; % beam drift evlocity in units of background thermal velocity
        B = 25;
        n = 0.06;
        R = 0.99;
        Te1 = 1600;
        Te2 = 60;
        Ti = 1600;
        n1= n*(1-R);
        n2= n*R;        
    case 2 % Daniels slow ESWs, lower beam temperature
        B = 25;
        vb_in_eV = 200;        
        n = 0.15;        
        n1 = 0.12;
        n2 = 0.03;
        R = n2/n;
        Te1 = 3500; 
        Te2  = 20; 
        Ti = 5000;
        S = sqrt(vb_in_eV/Te1);
    case 3 % Daniels fast ESWs
        B = 25;
        vb_in_eV = 1200;
        S = 0.2;
        
        n1 = 0.01;
        n2 = 0.003;
        n = n1+n2;        
        R = n2/n;
        Te1 = 4000; 
        Te2  = 250; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);
    case 4 % Daniels slow ESWs
        B = 25;
        vb_in_eV = 200;        
        n = 0.15;        
        n1 = 0.12;
        n2 = 0.03;
        R = n2/n;
        Te1 = 3500; 
        Te2  = 150; 
        Ti = 5000;
        S = sqrt(vb_in_eV/Te1);
        S = cn_eV2v(vb_in_eV,'eV')/cn_eV2v(Te1,'eV');
    case 42 % Daniels slow ESWs, lower beam temp
        B = 25;
        vb_in_eV = 300;        
        n = 0.15;        
        n1 = 0.12;
        n2 = 0.03;
        R = n2/n;
        Te1 = 3500; 
        Te2  = 30; 
        Ti = 5000;
        S = sqrt(vb_in_eV/Te1);
        S = cn_eV2v(vb_in_eV,'eV')/cn_eV2v(Te1,'eV');
    case 5 % Daniels fast ESWs, lower beam temp
        B = 25;
        vb_in_eV = 1200;
        S = 0.2;
        n = 0.013;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 100; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);   
    case 6 % Daniels fast ESWs, even lower beam temp
        B = 25;
        vb_in_eV = 1200;
        S = 0.2;
        n = 0.013;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);   
    case 7 % Daniels fast ESWs, even lower beam temp
        B = 25;
        vb_in_eV = 2000;
        S = 0.2;
        n = 0.1;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);  
    case 8 % Daniels fast ESWs, even lower beam temp and lower density
        B = 25;
        vb_in_eV = 3000;
        S = 0.2;               
        n1 = 0.01;
        n2 = 0.001;
        n = n1+n2; 
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1); 
    case 9 % Daniels fast ESWs, even lower beam temp and lower density
        B = 25;
        vb_in_eV = 3000;
        S = 0.2;               
        n1 = 0.009;
        n2 = 0.001;
        n = n1+n2; 
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);    
    case 10
        B = 25;
        vb_in_eV = 3000;
        S = 0.5;               
        n1 = 0.09;
        n2 = 0.01;
        n = n1+n2; 
        R = n2/n;
        Te1 = 1000; 
        Te2  = 600; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);    
    case 11 % Test case with small R
        B = 25;
        vb_in_eV = 3000;
        S = 0.5;               
        n=1;
        R=0.08;
        n1 = n*(1-R);
        n2 = n*R;
        %n1 = 0.09;
        %n2 = 0.01;
        %n = n1+n2; 
        %R = n2/n;
        Te1 = 500; 
        Te2  = Te1/25; 
        Ti = 500;
        S = 0.8;  
    case 12 % Daniels MP event
        %%
        B = 25;
        vb_in_eV = 90;
        S = 0.5;               
        n=1;
        R=0.08;
        n1 = n*(1-R);
        n2 = n*R;
        n1 = 1.8;
        n2 = 0.6;
        n = n1+n2; 
        R = n2/n;
        Te1 = 110; 
        Te2  = 17; 
        Ti = 300;
        S = 0.8;  
        S = sqrt(vb_in_eV/Te1)
    case 13 % Daniels MP event, tweaked beam
        %%
        B = 25;
        vb_in_eV = 70;
        S = 0.5;               
        n=1;
        R=0.08;
        n1 = n*(1-R);
        n2 = n*R;
        n1 = 1.8;
        n2 = 0.4;
        n = n1+n2; 
        R = n2/n;
        Te1 = 310; 
        Te2  = 11; 
        Ti = 310;
        S = 0.8;  
        S = sqrt(vb_in_eV/Te1)
    case 14 % Daniels MP event, tweaked beam
       %%
        B = 25;
        vb_in_eV = 90;
        S = 0.5;               
        n=1;
        R=0.08;
        n1 = n*(1-R);
        n2 = n*R;
        n1 = 1.8;
        n2 = 0.4;
        n = n1+n2; 
        R = n2/n;
        Te1 = 110; 
        Te2  = 6; 
        Ti = 300;
        S = 0.8;  
        S = sqrt(vb_in_eV/Te1)
    case 15 % Daniels new MP event
       %%
        B = 25;
        vb_in_eV = 80;
        n1 = 0.2+0.1;
        n2 = 0.3;
        n = n1+n2; 
        R = n2/n;
        Te1 = 350; 
        Te2  = 80; 
        Ti = 800;        
        S = sqrt(vb_in_eV/Te1)        
    case 16 % Daniels new MP event, lower beam temp
       %%
        B = 25;
        vb_in_eV = 80;
        n1 = 0.2+0.1;
        n2 = 0.3;
        n = n1+n2; 
        R = n2/n
        Te1 = 350; 
        Te2  = 14; 
        Ti = 800;        
        S = sqrt(vb_in_eV/Te1)
    case 162 % Daniels new MP event, lower beam temp
       %%
        B = 25;
        vb_in_eV = 80;
        n1 = 0.2+0.1;
        n2 = 0.3;
        n = n1+n2; 
        R = n2/n
        Te1 = 350; 
        Te2  = 14; 
        Ti = 300;        
        S = 0.48; 
    case 17 % compare to disprelESW_cn solver
       %%
        B = 25;       
        ni = 1; 
        R = 0.6;
        n1 = ni*(1-R);
        n2 = ni*R;
        
        Te1 = 300; 
        Te2  = 12; 
        Ti = 300;        
        S = 0.5;
    case 100 % Buneman marginal close to kLd=1.9
       %%
       % R=0.99; vb_in_eV = 63;
       % R=0.30; vb_in_eV = 46;
        B = 25;
        vb_in_eV = 46; 
        n=1;
        R=0.3;
        n1 = n*(1-R); 
        n2 = n*R;
        Te1 = 300; 
        Te2  = 12; 
        Ti = 300;        
        S = sqrt(vb_in_eV/Te1)  
        %SB = 1/sqrt(Te1/Te2)
    case 2540 % 
       %%
       % R=0.99; vb_in_eV = 63;
       % R=0.30; vb_in_eV = 46;
        B = 25;
        n=1;
        R=0.25;
        n1 = n*(1-R); 
        n2 = n*R;
        Te1 = 300; 
        Te2  = 12; 
        Ti = 300;        
        S = 0.4 ;
        %SB = 1/sqrt(Te1/Te2)
    case 2545 % 
       %%
       % R=0.99; vb_in_eV = 63;
       % R=0.30; vb_in_eV = 46;
        B = 25;
        n=1;
        R=0.25;
        n1 = n*(1-R); 
        n2 = n*R;
        Te1 = 300; 
        Te2  = 12; 
        Ti = 300;        
        S = 0.5 ;
        %SB = 1/sqrt(Te1/Te2)
    case 111 % fast esw MT
       %%
       if 0 % observed distribution
        B = 25;
        n1 = 0.01;
        n2 = 0.001;
        n = n1+n2;
        R = n2/n;         
        Te1 = 4000; 
        Te2  = 250; 
        Ti = 2000;  
        vb_in_eV = 1500; 
        S = 0.61;
        S = sqrt(vb_in_eV/Te1);
        vde1 = cn_eV2v(300,'eV')*1e3;
        %SB = 1/sqrt(Te1/Te2)
       else % tweaked
        B = 25;
        n1 = 0.01;
        n2 = 0.001;
        n = n1+n2;
        R = n2/n;         
        Te1 = 4000; 
        Te2  = 250; 
        Ti = 2000;  
        vb_in_eV = 1500; 
        S = 0.61;
        S = sqrt(vb_in_eV/Te1);
        vde1 = cn_eV2v(500,'eV')*1e3;
        %SB = 1/sqrt(Te1/Te2)
       end
end
%%

% Ions
omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
vdi = 0;
% Background electrons
omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
vthe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Vte'); % m/s 
%vde1 = 0; % m/s
% Beam electrons
omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
vthe2 = irf_plasma_calc(B,n,no,Te2,Ti,'Vte'); % m/s 
vde2 = S*vthe1;
% Other
lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m
omega_bune = omega_pe1^(1/3)*omega_pi^(2/3)*16^(-1/3);
omega_pe = irf_plasma_calc(B,n,no,Te2,Ti,'Fpe')*2*pi; % rad/s

normal_frequency = omega_pe1;
normal_length = lamD;
normal_frequency_string = '\omega_{pe}';
normal_length_string = '\lambda_D';


%vde1= 0*(-vde2);
fv = [omega_pi,omega_pe1,omega_pe2,vthi,vthe1,vthe2,vdi,vde1,vde2,lamD];

% set up figure for drawing results
figure(3)    
for nPlot = 1:4   
    h(nPlot) = subplot(1,4,nPlot); hold(h(nPlot),'on');        
    xlabel(h(nPlot),['k' normal_length_string],'Fontsize',14);
    ylabel(h(nPlot),['\omega/'  normal_frequency_string],'Fontsize',14);    
    set(h(nPlot),'XGrid','on','YGrid','on','Fontsize',14);
end
ylabel(h(3),'v_{ph} [km/s]')
title(h(1),'\Re (\omega)')
title(h(2),'\Im (\omega)')
title(h(3),'v_{ph}')
xlabel(h(4),'v [10^3 km/s]')
title(h(4),'Input distribution')
v_lim = 2*6*cn_eV2v(Te2,'eV');
v_psd = -v_lim:100:v_lim;
Units.eV = 1.6022e-19;
Units.me = 9.109399999999e-31;
eV_psd = Units.me*v_psd.^2*10^6/2/Units.eV;
e_psd1 = cn.maxwellian(v_psd,Te1,n1,vde1/1000,'e');
e_psd2 = cn.maxwellian(v_psd,Te2,n2,vde2/1000,'e');
i_psd = cn.maxwellian(v_psd,Ti,n1+n2,vdi/1000,'p');
e_psd = e_psd1+e_psd2;

plot(h(4),v_psd,e_psd,v_psd,e_psd1,v_psd,e_psd2,v_psd,i_psd/i_psd*max(e_psd),'k--')
%set(h(4),'yscale','log','xscale','log','xlim',[1e1 2e4])
set(h(4),'xlim',v_lim*[-1 1])
axes(h(4))
distrStr = {['R=' num2str(R,'%.2f')],...
            ['Ti=' num2str(Ti,'%.0f')],...
            ['Te1=' num2str(Te1,'%.0f')],...
            ['Te2=' num2str(Te2,'%.0f')],...
            ['S=' num2str(S,'%.2f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')

%%

ks = 1:0.05:1.8;
ks = 0.2:0.05:1.8;
ks = 0.01:0.05:1;
ks = 0.1:0.05:1.5;
ks = 0.2:0.05:1;
ks = 1.04:0.01:1.1;
ks = 0.1:0.05:1.3;
%ks = 0.481:0.001:0.49;
%ks = 1.05:0.05:1.4;
rights = zeros(1,numel(ks));
nTries = 6;
wimax_store = nan(1,numel(ks));
wrmax_store = nan(1,numel(ks));
for na=1:numel(ks);    
    a = ks(na);
    wimax = 0;
    % For each k*lD, make 20 random guesses within specified ranges
    for iii=1:nTries
        if iii<fix(nTries/2+1) % (nTries+1)
            x_real=0.1*randn*normal_frequency;             
            x_imag=0.05*randn*normal_frequency; % x_imag=2*rand*omega_pe;  
        else % around normal_frequency              
            x_real=(1*randn*normal_frequency);            
            x_imag=0.2*randn*normal_frequency; % x_imag=2*rand*omega_pe;
        end
        %af = @(temp) fff_two_stream(temp,a);
        %af = @(temp) fff_two_stream_test(temp,a);
        af = @(temp) myfun_20070831_EH(temp,a,fv);
        %[x,FVAL,EXITFLAG] = lsqnonlin(af,[x_real x_imag],0+0i,10000+0i,optimoptions('lsqnonlin','MaxIter',2000,'TolFun',TolFun,'TolX',1e-16,'MaxFunEvals',5000,'display','off'));
        [x,FVAL,EXITFLAG] = fsolve(af,[x_real x_imag],optimoptions('fsolve','MaxIter',1000,'TolFun',1e-10,'TolX',1e-13,'MaxFunEvals',5000,'display','off'));
        [x_real, x_imag];
        if(EXITFLAG==1) % Found something, now plot            
            plot(h(1),real(a),x(1)/normal_frequency,'k.','Markersize',3);             
            plot(h(1),real(a),x(2)/normal_frequency,'k.','Markersize',3); hold on;
            plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'K.','Markersize',3); hold on; 
            % only plot for reasonably omega_i > -0.6
            if(x(2)/normal_frequency>-0.6)
                if ( abs(x(1)/normal_frequency)<0.8)                    
                    plot(h(1),real(a),x(1)/normal_frequency,'r+','Markersize',6); hold on; 
                    plot(h(2),real(a),x(2)/normal_frequency,'r+','Markersize',6); hold on; 
                    plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'R+','Markersize',6); hold on; 
                elseif 0;( abs(x(1)/normal_frequency)>0.8) % before it was 0.1                    
                    plot(h(1),real(a),x(1)/normal_frequency,'bd','Markersize',6); hold on;                    
                    plot(h(2),real(a),x(2)/normal_frequency,'bd','Markersize',6); hold on;    
                end
                if (x(2)/normal_frequency>0.0001)  % positive growth rate ?                  
                    rights(na)=rights(na)+1;
                    if x(2) > wimax
                        wimax_store(na) = x(2);
                        wrmax_store(na) = x(1);
                    end
                    plot(h(1),real(a),x(1)/normal_frequency,'ch','Markersize',6); hold on;                    
                    plot(h(2),real(a),x(2)/normal_frequency,'ch','Markersize',6); hold on; 
                    plot(h(3),real(a),x(1)/real(a)*normal_length/1000,'ch','Markersize',6); hold on; 
                    disp(['x_real=' num2str(x_real,'%.4f') ', x_imag=' num2str(x_imag,'%.4f') ' => x_real=' num2str(x(1),'%.4f') ', x_imag=' num2str(x(2),'%.4f')])
                end
            end                   
            drawnow               
        end
    end
    rights;
    % Some analytical solution?
    %plot(h(1),real(a),sqrt(omega_pe*omega_pe + 3*(a/lamD)*(a/lamD)*vthe*vthe/2)/omega_pe  ,'r.','Markersize',3); hold on;
    %plot(h(1),real(a),(vSound*real(a)/2/lamD)/omega_pe,'g.','Markersize',3); hold on;   
end

if 0
    set(h(1),'ylim',[0 0.6],'xlim',[0 max(ks)]);
    set(h(2),'ylim',[-0.01 0.04],'xlim',[0 max(ks)])
end

if 1 % add phase velocities to distribution plot
    %%
    nvs = input('How many vph would you like to add?')
    hold(h(4),'on')
    %set(h(3),'ylim',1e3*[-30 30])
    for kk = 1:nvs
        vph_toplot = ginput(1);
        plot(h(4),vph_toplot(2)*[1 1],[0 max(e_psd)],'r--')
    end
    %set(h(4),'ylim',[0 max(e_psd)])
    hold(h(4),'off')
end
    