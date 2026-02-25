% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

% Distribution properties
ni = 1;
R = 0.25;
ne1 = ni*(1-R);
ne2 = ni*R;
Te1 = 300;
Te2 = 12;
Ti = 300;
S=0.7;

qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

vith = sqrt(2*qe*Ti/mi);
veth1 = sqrt(2*qe*Te1/me);
veth2 = sqrt(2*qe*Te2/me);
vd1 = veth1*0.0; %vd1 = sqrt(2*qe*300/me);
vd2 = veth1*S; %vd2 = sqrt(2*qe*1500/me);

wpe = sqrt((ne1+ne2)*qe^2/(me*eps)); % Hz
wpe1 = sqrt(ne1*qe^2/(me*eps)); % Hz
wpe2 = sqrt(ne2*qe^2/(me*eps)); % Hz
wpi = sqrt(ni*qe^2/(mi*eps)); % Hz
Ld = sqrt(veth1^2)/wpe/sqrt(2);


nk = 200;
kvec = linspace(1e-3,1,nk)/Ld;
dk = kvec(2)-kvec(1);
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));


%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
% higher vde2 -> higher vph, this should control the guessing range
% w/k = vph -> wrange >~vph*k

wp1 = wpe1;
wp2 = wpe2;
wp3 = wpi;
vt1 = veth1;
vt2 = veth2;
vt3 = vith;
vd1 = vd1;
vd2 = vd2;
vd3 = 0;

wapprox = kvec(1)*vd2;
wrvec = linspace(0,wapprox*0.5,1501);
wivec = linspace(0,wapprox*0.5,1401);
[wr,wi] = meshgrid(wrvec,wivec);
wc = wr+i*wi;

% Dispersion relation
Nfad = 50;
xi = @(w,k,vd,vt) (w-k*vd)/(k*vt);
psi = @(w,k,vd,vt,wp) 2*wp^2/(k^2*vt^2)*(1 + i*sqrt(pi)*xi(w,k,vd,vt).*faddeeva(xi(w,k,vd,vt),Nfad));
D = @(w,k) 1 + psi(w,k,vd1,vt1,wp1) + psi(w,k,vd2,vt2,wp2) + psi(w,k,vd3,vt3,wp3);

z = abs(D(wc,kvec(1))); %min(min(z));
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
plot(h(1),kvec(1)*Ld,wr_(1)/wpi,'o',kvec(1)*Ld,wi_(1)/wpi,'x'); hold(gca,'on')
    
frange = 1.5*wapprox;
for nn = 2:length(kvec)
  frange = vd2*dk*0.5;
  %linear search for zero point around linearly interpolated guess. frange
  %may need to be changed depending on mode. Large frange decreases accuracy
  nf = 401;
  wrv = linspace(wr_(nn-1)+gradr-frange,wr_(nn-1)+gradr+frange,nf);
  wiv = linspace(wi_(nn-1)+gradi-frange,wi_(nn-1)+gradi+frange,nf);
  %if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end
  wrv(wrv<0)=[];

  [wr,wi] = meshgrid(wrv,wiv);
  wc = wr+i*wi;

  z = abs(D(wc,kvec(nn)));
  min(min(z))
  [q r] = find(min(min(z))==z);
  wr_(nn) = wrv(r);
  wi_(nn) = wiv(q);

  wnorm = 1; wpi;
  plot(h(1),kvec(nn)*Ld,wr_(nn),'o',kvec(nn)*Ld,wi_(nn)/wnorm,'x')
  errorbar(h(1),kvec(nn)*Ld,wrv(fix(nf/2))/wnorm,frange/wnorm)
  errorbar(h(1),kvec(nn)*Ld,wiv(fix(nf/2))/wnorm,frange/wnorm)
  pause(0.01)
  %To do: Write better root finder.

  gradr = wr_(nn)-wr_(nn-1);
  gradi = wi_(nn)-wi_(nn-1);
end


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