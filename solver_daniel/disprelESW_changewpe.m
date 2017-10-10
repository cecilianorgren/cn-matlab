% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

% Distribution properties 
ne1 = 0.05e6;
ne2 = 0.55e6;
ni = 0.6e6;
Te1 = 350;
Te2 = 14;
Ti = 800;

R = 0.9;
S = 0.5;
ni = 1;
ne1 = ni*(1-R);
ne2 = ni*R;
Te1 = 300;
Te2 = 12;
Ti = 300;

qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps = 8.854e-12;
kb = 1.381e-23;

vith = sqrt(2*qe*Ti/mi);
veth1 = sqrt(2*qe*Te1/me);
veth2 = sqrt(2*qe*Te2/me);
vd1 = veth1*0.0;
vd2 = veth1*S;


wpe = sqrt((ne1+ne2)*1e6*qe^2/(me*eps)); % Hz
wpe1 = sqrt(ne1*1e6*qe^2/(me*eps));
wpe2 = sqrt(ne2*1e6*qe^2/(me*eps));
wpi = sqrt(ni*1e6*qe^2/(mi*eps));
Ld = sqrt(veth1^2)/wpe/sqrt(2); % m

dk = 2e-5;
kvec = [2e-5:dk:0.008];
wia = length(kvec);
wga = length(kvec);


%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
wrv = linspace(0,10,1001);
wiv = linspace(0,10,1001);
[wr,wi] = meshgrid(wrv,wiv);
wc = wr+i*wi;


    xie1 = @(omega)((omega-kvec(1)*vd1)/(kvec(1)*veth1));
    xie2 = @(omega)((omega-kvec(1)*vd2)/(kvec(1)*veth2));
    xii = @(omega)(omega/(kvec(1)*vith));
    iadisp = @(w) (1+2*wpe1^2/(kvec(1)^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w).*faddeeva(xie1(w),50))+...
                     2*wpe2^2/(kvec(1)^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w).*faddeeva(xie2(w),50))+...
                     2*wpi^2/(kvec(1)^2*vith^2)*(1+i*sqrt(pi)*xii(w).*faddeeva(xii(w),50)));
    z = abs(iadisp(wc));
    min(min(z));
    [q r] = find(min(min(z))==z);
    wia(1) = wrv(r);
    wga(1) = wiv(q);

    gradr = wia(1);
    gradi = wiv(1);

    h = subplot(1,1,1);
plot(h(1),kvec(1)*Ld,wia(1)/wpi,'o',kvec(1)*Ld,wga(1)/wpi,'x'); hold(gca,'on')
    
    
frange = 8;
for nn = 2:length(kvec)
%linear search for zero point around linearly interpolated guess. frange
%may need to be changed depending on mode. Large frange decreases accuracy
wrv = linspace(wia(nn-1)+gradr-frange,wia(nn-1)+gradr+frange,401);
wiv = linspace(wga(nn-1)+gradi-frange,wga(nn-1)+gradi+frange,401);
if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end

[wr,wi] = meshgrid(wrv,wiv);
wc = wr+i*wi;


    xie1 = @(omega)((omega-kvec(nn)*vd1)/(kvec(nn)*veth1));
    xie2 = @(omega)((omega-kvec(nn)*vd2)/(kvec(nn)*veth2));
    xii = @(omega)(omega/(kvec(nn)*vith));
    iadisp = @(w) (1+2*wpe1^2/(kvec(nn)^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w).*faddeeva(xie1(w),50))+...
                     2*wpe2^2/(kvec(nn)^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w).*faddeeva(xie2(w),50))+...
                     2*wpi^2/(kvec(nn)^2*vith^2)*(1 + i*sqrt(pi)*xii(w).*faddeeva(xii(w),50)));

    z = abs(iadisp(wc));
    min(min(z))
    [q r] = find(min(min(z))==z);
    wia(nn) = wrv(r);
    wga(nn) = wiv(q);
    
    plot(h(1),kvec(nn)*Ld,wia(nn)/wpi,'o',kvec(nn)*Ld,wga(nn)/wpi,'x')
    errorbar(h(1),kvec(nn)*Ld,wrv(fix(nf/2))/wpi,frange/wpi)
    errorbar(h(1),kvec(nn)*Ld,wiv(fix(nf/2))/wpi,frange/wpi)
    pause(0.01)
    
    %To do: Write better root finder.
    
    gradr = wia(nn)-wia(nn-1);
    gradi = wga(nn)-wga(nn-1);
end


dispcurv = struct;
dispcurv.kvec = kvec;
dispcurv.wr = wia;
dispcurv.wi = wga; 
dispcurv.wpi = wpi; 
dispcurv.wpe = wpe;

if 0,
save('modbun.mat','dispcurv');
end

pos=find(max(wga)==wga)
kmax = kvec(pos(1))
vph = wia(pos(1))/kvec(pos(1))
vphnorm = vph(1)/veth1
maxfreq = wia(pos(1))/wpi
maxgamma = wga(pos(1))/wpi

plot(kvec,wia/wpi,kvec,wga/wpi); xlabel('k [m^{-1}]')
plot(kvec*Ld,wia/wpi,kvec*Ld,wga/wpi); xlabel('k\lambda_{De}')

ylabel('\omega_{max}/\omega_{pi}, \gamma_{max}/\omega_{pi}')
%find(min(abs(iadisp(wc))))