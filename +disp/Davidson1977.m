% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

units = irf_units;

% Distribution properties
n = 4e6;
ne = n;
ni = n;
Te = 150; Te_K = Te*units.eV/units.kB;
Ti = 1000; Ti_K = Ti*units.eV/units.kB;
B = 40e-9;

qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

vti = sqrt(2*qe*Ti/mi);
vte = sqrt(2*qe*Te/me);
wpe = sqrt(ne*qe^2/(me*eps)); % Hz
wpi = sqrt(ni*qe^2/(mi*eps)); % Hz
wce = qe*B/me;
wci = qe*B/mi;
roe = vte/wce;
roi = vti/wci;
Ld = sqrt(vte^2)/wpe/sqrt(2)*1e-3;

wlh = sqrt(wce*wci);
% Gradient length scales
Ln = 1/(2e-5);5*roe;10e3;
LB = 10e5;
LT = 10e5;

nk = 80;
kvec = 2*linspace(1e-3,1,nk)/roe;
dk = kvec(2)-kvec(1);
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));

vdi = 0;


% Dispersion relation
Nfad = 50;
xi = @(w,k,vd,vt) (w-k*vd)/(k*vt);
%psi_i = @(w,k,vd,vt,wp) 2*wp^2/(k^2*vt^2)*(1 + i*sqrt(pi)*xi(w,k,vd,vt).*faddeeva(xi(w,k,vd,vt),Nfad));
b = @(k,vt,wc) k.^2*vt^2/2/wc^2;
% v_delta = @(k,vt,wc) ...
%     - 0.5*vt.^2/wce*besseli(0,b(k,vt,wc)).*exp(-besseli(0,b(k,vt,wc))).*(...
%         1/Ln-...
%         1/LB*(1-b(k,vt,wc).*(1-besseli(1,b(k,vt,wc))./besseli(0,b(k,vt,wc))))-...
%         1/LT*b(k,vt,wc).*(1-besseli(1,b(k,vt,wc))./besseli(0,b(k,vt,wc))));
v_delta = @(k) ...
    + 0.5*vte.^2/wce*besseli(0,b(k,vte,wce)).*exp(-b(k,vte,wce)).*(...
        1/Ln);
vE = 600e3;

%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
% higher vde2 -> higher vph, this should control the guessing range
% w/k = vph -> wrange >~vph*k

vguess = 100e3;vE;sqrt(wci*wce)*roe;
wapprox = wlh;kvec(1)*vguess;
wrvec = linspace(-1,1,1501);
wivec = linspace(-1,1,1401);
[wr,wi] = meshgrid(wrvec,wivec);
wmat = wr+i*wi;


vdi = 0;

psi_i = @(w,k,vd,vt,wp,wc) 2*wp^2./(k.^2)/(vti^2).*(1+w./k/vti.*i*sqrt(pi).*faddeeva(w./k/vti,Nfad));
psi_e = @(w,k,vt,wp,wc) (wp^2/wc^2).*(1-besseli(0,b(k,vt,wc)).*exp(-b(k,vt,wc)))./b(k,vt,wc) + ...
                           2*wp^2./k.^2/vt^2.*k.*v_delta(k)./(w-k*vE);

D = @(w,k) 1 + psi_i(w,k,vdi,vti,wpi) + psi_e(w,k,vte,wpe,wce);

z = abs(D(wmat,kvec(1))); %min(min(z));
%pcolor(abs(z)); shading flat
[q r] = find(min(min(z))==z);
wr_(1) = wrvec(r);
wi_(1) = wivec(q);
 
gradr = wr_(1);
gradi = wi_(1);

if 0 % check grid    
    %%
    subplot(3,1,1)
    pcolor(wr,wi,real(D(wc)))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,2)
    pcolor(wr,wi,imag(D(wc)))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,3)
    pcolor(wr,wi,abs(D(wc)))
    xlabel('w_r')
    ylabel('w_i')
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[0 1])
end
  
h = subplot(1,1,1);
plot(h(1),kvec(1)*roe,wr_(1)/wlh,'o',kvec(1)*roe,wi_(1)/wlh,'x'); hold(gca,'on')
    
frange = 1.5*wapprox;
for nn = 2:length(kvec)
  frange = vguess*dk*2;
  %linear search for zero point around linearly interpolated guess. frange
  %may need to be changed depending on mode. Large frange decreases accuracy
  nf = 201;
  wrv = linspace(wr_(nn-1)+gradr-frange,wr_(nn-1)+gradr+frange,nf);
  wiv = linspace(wi_(nn-1)+gradi-frange,wi_(nn-1)+gradi+frange,nf);
  %if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end
  %wrv(wrv<0)=[];

  [wr,wi] = meshgrid(wrv,wiv);
  wc = wr+i*wi;

  z = abs(D(wc,kvec(nn)));
  ze = abs(psi_e(wc,kvec(nn),vte,wpe,wce));
  zi = abs(psi_i(wc,kvec(nn),vdi,vti,wpi,wci));
  %z=ze;
  min(min(z))
  [q r] = find(min(min(z))==z);
  wr_(nn) = wrv(r);
  wi_(nn) = wiv(q);

  wnorm = wlh; wpi;
  if 1
    plot(h(1),kvec(nn)*roe,wr_(nn)/wnorm,'o',kvec(nn)*roe,wi_(nn)/wnorm,'x')
    %errorbar(h(1),kvec(nn)*roe,wrv(fix(nf/2))/wnorm,frange/wnorm)
    %errorbar(h(1),kvec(nn)*roe,wiv(fix(nf/2))/wnorm,frange/wnorm)
  else
    pcolor(abs(z)); shading flat;
    colorbar
    set(gca,'clim',[0 1e2])
  end
  pause(0.1)
  %To do: Write better root finder.

  gradr = wr_(nn)-wr_(nn-1);
  gradi = wi_(nn)-wi_(nn-1);
end

%%
dispcurv = struct;
dispcurv.kvec = kvec;
dispcurv.wr = wr_;
dispcurv.wi = wi_; 
dispcurv.wpi = wpi; 
dispcurv.wpe = wpe;
dispcurv.R = 0.25;
dispcurv.ne1 = ni*(1-R);
dispcurv.ne2 = ni*R;
dispcurv.Te1 = 300;
dispcurv.Te2 = 12;
dispcurv.Ti = 300;
dispcurv.S = 0.5;
dispcurv.Ld = Ld;
dispcurv.veth1 = veth1;
dispcurv.veth2 = veth2;
dispcurv.vith = vith;

if 0,
save('modbun.mat','dispcurv');
end

pos=find(max(wga)==wga);
kmax = kvec(pos(1));
vph = wia(pos(1))/kvec(pos(1));
vphnorm = vph(1)/veth1;
maxfreq = wia(pos(1))/wpi;
maxgamma = wga(pos(1))/wpi;

if 0
%%

plot(kvec*Ld,wr_/wpi,kvec*Ld,wi_/wpi)
plot(kvec*Ld,wr_/wpe,kvec*Ld,wi_/wpe)
xlabel('k\lambda_D')
ylabel('\omega_{max}/\omega_{pi}, \gamma_{max}/\omega_{pi}')
distrStr = {['R=' num2str(ne2/ni,'%.2f')],...
            ['Ti=' num2str(Ti,'%.0f')],...
            ['Te1=' num2str(Te1,'%.0f')],...
            ['Te2=' num2str(Te2,'%.0f')],...
            ['S=' num2str(vd2/veth1,'%.2f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')
end