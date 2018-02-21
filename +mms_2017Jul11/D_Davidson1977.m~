function D = D_Davidson1977(w,k,vte,wce,wpe,vti,wci,wpi,vE,LB,Ln,LT,option,printD)
c = 299792458;

vB = -1/LB*vte^2/2/wce;
vn = -1/Ln*vte^2/2/wce;
vT = -1/LT*vte^2/2/wce;

Nfad = 50;
nv = 100;
vvec = linspace(1,5*vte,nv);
ww = repmat(w,1,1,nv);
vmat = permute(repmat(vvec,size(w,1),1,size(w,2)),[1 3 2]);

%psi_i = @(w,k,vt,wp,wc) 2*wp^2./(k.^2)/(vti^2).*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad));

if ~exist('option','var') || isempty(option), option = 'full'; end
%disp(option)
switch option
  case 'colde' % cold electrons, Davidson1977 (48) with (50) substitutions    
    % Ln<0
    % LT<0
    % LB>0
    D = + 1 ...
        + 2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad)) ... % ions      
        + wpe^2/wce^2 ...
        - wpe^2/k/wce./(w-k*vE)*(1/Ln-1/LB) ... 
        + 2*wpe^4/c^2/k^2/wce^2;
  case 'colde_noEM' % cold electrons, Davidson1977 (48) with (50) substitutions
    % the imaginary terms, i.e. everything containing w, are ok
    % d77 just approximates I0(b)*exp(-b)~0, this makes just a small
    % difference
    D = + 1 ...
        + 2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad)) ... % ions      
        + wpe^2/wce^2 ...
        - wpe^2/k/wce./(w-k*vE)*(1/Ln-1/LB) ... 
        + 0*2*wpe^4/c^2/k^2/wce^2 ...
        ;              
  case 'reduced' % Reduced D from Davidson1977, to debug and compare to 75
    % besselj: Bessel function of the first kind
    J0p = @(x) (0*besselj(0, x))./x - besselj(0 + 1, x); % x is mu = k*v/wce
    J0 = @(x) besselj(0, x); % x is mu = k*v/wce

    A = @(w,k,v) w - k*vE - k*vn - k*vT*(1-v.^2/vte^2);
    %A = @(w,k,v) w - k*vE;

    O1_ = @(w,k,v) 2/vte^4*v.^3.*J0p(k*v/wce).^2          .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
    O2_ = @(w,k,v) 2/vte^3*v.^2.*J0(k*v/wce).*J0p(k*v/wce).*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
    O3_ = @(w,k,v) 2/vte^2*v.^1.*J0(k*v/wce).^2           .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);


    b = @(wce,k,vt) k.^2*vt^2/2./wce.^2;
    O3_comp_ = @(wce,k,vt) besseli(0,b(wce,k,vt)).*exp(-b(wce,k,vt));
    O3_comp = O3_comp_(wce,k,vte);

    O1 = sum(O1_(ww,k,vmat),3)*(vvec(2)-vvec(1));
    O2 = sum(O2_(ww,k,vmat),3)*(vvec(2)-vvec(1));
    O3 = sum(O3_(ww,k,vmat),3)*(vvec(2)-vvec(1));
    %O2 = O3*0;
    %O1 = O3*0;
    %O3 = O3_comp;

    %O1 = integral(@(v)O1_(w(100,100),k,v),0,Inf)
    %O1 = trapz(O1_(w(100,100),k,vvec),vvec)
    %O2 = trapz(O2_(w,k,vvec),vvec);
    %O3 = trapz(O3_(w,k,vvec),vvec);

    D = + (2*wpe^2/c^2/k^2)*(2*wpe^2/k^2/vte^2)*O2.^2 ...
        + (1+2*wpe^2/c^2/k^2*O1).*(1+2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad)) ...
            + 2*wpe^2/k^2/vte^2*(1-O3));
  case 'full' % Full D from Davidson1977
    % besselj: Bessel function of the first kind
    J0p = @(x) (0*besselj(0, x))./x - besselj(0 + 1, x); % x is mu = k*v/wce
    J0 = @(x) besselj(0, x); % x is mu = k*v/wce

    A = @(w,k,v) w - k*vE - k*vn - k*vT*(1-v.^2/vte^2);
    %A = @(w,k,v) w - k*vE;

    O1_ = @(w,k,v) 2/vte^4*v.^3.*J0p(k*v/wce).^2          .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
    O2_ = @(w,k,v) 2/vte^3*v.^2.*J0(k*v/wce).*J0p(k*v/wce).*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);
    O3_ = @(w,k,v) 2/vte^2*v.^1.*J0(k*v/wce).^2           .*exp(-v.^2/vte^2)./(w-k*vE-k*vB).*A(w,k,v);

    %b = @(wce,k,vt) k.^2*vt^2/2./wce.^2;
    %O3_comp_ = @(wce,k,vt) besseli(0,b(wce,k,vt)).*exp(-b(wce,k,vt));
    %O3_comp = O3_comp_(wce,k,vte);
    
    % Integration over vperp 'zero to infty'
    O1 = sum(O1_(ww,k,vmat),3)*(vvec(2)-vvec(1));
    O2 = sum(O2_(ww,k,vmat),3)*(vvec(2)-vvec(1)); % O2 only appears as squared
    O3 = sum(O3_(ww,k,vmat),3)*(vvec(2)-vvec(1));
    
    %O1 = integral(@(v)O1_(w(100,100),k,v),0,Inf)
    %O1 = trapz(O1_(w(100,100),k,vvec),vvec)
    %O2 = trapz(O2_(w,k,vvec),vvec);
    %O3 = trapz(O3_(w,k,vvec),vvec);

    D = ...
        + (1+2*wpe^2/c^2/k^2*O1).* ...
          (1+2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad)) + 2*wpe^2/k^2/vte^2*(1-O3)) ...
        + (2*wpe^2/c^2/k^2)*(2*wpe^2/k^2/vte^2)*O2.^2 ...    
        ;
     if 0%printD
     fprintf('Di = (1+%g+%g)*%g+%g    ',imag(2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad))),imag(2*wpe^2/k^2/vte^2*(1-O3)),imag((1+2*wpe^2/c^2/k^2*O1)),imag((2*wpe^2/c^2/k^2)*(2*wpe^2/k^2/vte^2)*O2.^2)) 
     fprintf('Dr = (1+%g+%g)*%g+%g\n',real(2*wpi^2/k^2/vti^2*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad))),real(2*wpe^2/k^2/vte^2*(1-O3)),real((1+2*wpe^2/c^2/k^2*O1)),real((2*wpe^2/c^2/k^2)*(2*wpe^2/k^2/vte^2)*O2.^2)) 
     end
     % compare to colde approximation
     O1_approx = 0;
     O2_approx = k*vte/2/wce;
     O3_approx = 1 - k^2*vte^2/2/wce^2 + k*vte^2/2/wce./(w-k*vE)*(1/Ln-1/LB);               
end