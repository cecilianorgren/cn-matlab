k1=0.08; k2=2; nk=100; k=linspace(k1,k2,nk);

zeta=@(x)faddeeva(x)*1i*sqrt(pi);

f0 = @(o,k,Ti,Te,Tb,vb,vteb,vte,vi,me,mi,R)(...
        k^2 + ...
        (Ti/Te)*(me/mi)^2*(1+o/k*sqrt(me*Ti/mi/Te/2).*zeta(o/k*sqrt(me*Ti/mi/Te/2))) + ...
        (1-R)*(1+o/k*sqrt(1/2).*zeta(o/k*sqrt(1/2))) + ...
        R*(Tb/Te)*(1+(o/k*sqrt(Tb/Te/2)-vb/vteb).*zeta((o/k*sqrt(Tb/Te/2)-vb/vteb))));
        
% Model 
units=irf_units;
R = 1;  % beam ratio
Ti= 2000; % eV
Te= 1600; % eV
Tb= 60;   % eV
vteb=4600e3; % m/s
vb= 5.5 * vteb; % m/s
vte =24000e3; % m/s
vi =  760e3; % m/s
mi = units.mp;
me = units.me;
    
    
    

f = @(o,k) f0(o,k,Ti,Te,Tb,vb,vteb,vte,vi,me,mi,R);
w=[];
for kk=k
    options=optimset('Display','off'); 
    om=fsolve(f,1e2+0.1i,options,kk)*sqrt(2)*kk; 
    w=[w,om];
end

wi = imag(w);
wr = real(w);
h(1)=subplot(1,2,1); plot(h(1),k,wr,'*r')
set(h(1),'xlim',k([1 end]),'yscale','log')
ylabel(h(1),'\omega_r/\omega_{pe}')
xlabel(h(1),'x\lambda_{de}')
%legend(h(1),'\omega_r','\omega_i')
h(2)=subplot(1,2,2); plot(h(2),k,wi,'*b')
set(h(2),'xlim',k([1 end]))
ylabel(h(2),'\omega_i/\omega_{pe}')
xlabel(h(2),'x\lambda_{de}')

    