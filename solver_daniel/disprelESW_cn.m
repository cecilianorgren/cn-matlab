% 3-component dispersion equation solver - D. B. Graham
% When running if ans does not go to low values [below 1e-(3--4)] the
% solution is not valid.

% Distribution properties 
%ne1 = 0.01e6;
%ne2 = 0.003e6;
%ni = ne1+ne2;0.013e6;
ni = 1;
R = 0.25;
ne1 = ni*(1-R);
ne2 = ni*R;
Te1 = 300;
Te2 = 12;
Ti = 300;
S=0.7;

if 0
ne1 = 0.3e6;
ne2 = 0.3e6;
ni = 0.6e6;
Te1 = 350;
Te2 = 14;
Ti = 800;
end
if 0
    % Magnetoail, fast holes
    B = 25;
    vb_in_eV = 1200;                
    ne1 = 0.01;
    ne2 = 0.003;
    ni = ne1+ne2;        
    R = ne2/ni;
    Te1 = 4000; 
    Te2  = 250; 
    Ti = 2000;
    S = sqrt(vb_in_eV/Te1);
end
if 0
    % Magnetoail, slow holes
    B = 25;
    vb_in_eV = 200;                
    ne1 = 0.12;
    ne2 = 0.03;
    ni = n1+n2;        
    R = n2/n;
    Te1 = 3500; 
    Te2  = 150; 
    Ti = 5000;
    S = sqrt(vb_in_eV/Te1);
end
if 0
ne1 = 0.01e6;
ne2 = 0.001e6;
ni = 0.011e6;
Te1 = 4000;
Te2 = 30;
Ti = 2000;
end

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




dk = 2e-5;
kvec1 = [2e-1:dk:0.008];
nk = 200;
kvec2 = linspace(1e-3,1,nk)/Ld;
kvec = kvec2;
dk = kvec(2)-kvec(1);
wia = nan(1,length(kvec));
wga = nan(1,length(kvec));


%Initial range of complex frequencies at low k. Linspace range may need to
%be varied depending on mode.
% higher vde2 -> higher vph, this should control the guessing range
% w/k = vph -> wrange >~vph*k

wapprox = kvec(1)*vd2;
wrv = linspace(0,wapprox*0.5,1501);
wiv = linspace(0,wapprox*0.5,1401);
%wiv(wiv<0)=[];
%wrv(wrv<0)=[];
[wr,wi] = meshgrid(wrv,wiv);
wc = wr+i*wi;


    xie1 = @(omega)((omega-kvec(1)*vd1)/(kvec(1)*veth1));
    xie2 = @(omega)((omega-kvec(1)*vd2)/(kvec(1)*veth2));
    xii = @(omega)(omega/(kvec(1)*vith));
    iadisp = @(w) (1+2*wpe1^2/(kvec(1)^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w).*faddeeva(xie1(w),50))+...
                     2*wpe2^2/(kvec(1)^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w).*faddeeva(xie2(w),50))+...
                     2*wpi^2/(kvec(1)^2*vith^2)*(1 + i*sqrt(pi)*xii(w).*faddeeva(xii(w),50)));
    z = abs(iadisp(wc));
    min(min(z));
    [q r] = find(min(min(z))==z);
    wia(1) = wrv(r);
    wga(1) = wiv(q);

    gradr = wia(1);
    gradi = wiv(1);
%
if 0 % check grid    
    %%
    subplot(3,1,1)
    pcolor(wr,wi,real(iadisp(wc)))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,2)
    pcolor(wr,wi,imag(iadisp(wc)))
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[-1 1])
    subplot(3,1,3)
    pcolor(wr,wi,abs(iadisp(wc)))
    xlabel('w_r')
    ylabel('w_i')
    shading('flat')
    colorbar
    %set(gca,'clim',1e1*[0 1])
end
%  
h = subplot(1,1,1);
plot(h(1),kvec(1)*Ld,wia(1)/wpi,'o',kvec(1)*Ld,wga(1)/wpi,'x'); hold(gca,'on')
    
xie1 = @(omega,k)((omega-k*vd1)/(k*veth1));
xie2 = @(omega,k)((omega-k*vd2)/(k*veth2));
xii = @(omega,k)(omega/(k*vith));
iadisp = @(w,k) (1+2*wpe1^2/(k^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w,k).*faddeeva(xie1(w,k),50))+...
                 2*wpe2^2/(k^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w,k).*faddeeva(xie2(w,k),50))+...
                 2*wpi^2/(k^2*vith^2)*(1 + i*sqrt(pi)*xii(w,k).*faddeeva(xii(w,k),50)));

             
%
frange = 5;
frange = 1.5*wapprox;
%2e-3*wpe;
for nn = 2:length(kvec)
    frange = vd2*dk*0.5;
%linear search for zero point around linearly interpolated guess. frange
%may need to be changed depending on mode. Large frange decreases accuracy
nf = 401;
wrv = linspace(wia(nn-1)+gradr-frange,wia(nn-1)+gradr+frange,nf);
wiv = linspace(wga(nn-1)+gradi-frange,wga(nn-1)+gradi+frange,nf);
%if kvec(nn)*Ld<0.2; wiv(wiv<0)=[]; end
wrv(wrv<0)=[];

[wr,wi] = meshgrid(wrv,wiv);

wc = wr+i*wi;
%wc(wi<0)=[];

if 0
    xie1 = @(omega)((omega-kvec(nn)*vd1)/(kvec(nn)*veth1));
    xie2 = @(omega)((omega-kvec(nn)*vd2)/(kvec(nn)*veth2));
    xii = @(omega)(omega/(kvec(nn)*vith));
    iadisp = @(w) (1+2*wpe1^2/(kvec(nn)^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w).*faddeeva(xie1(w),50))+...
                     2*wpe2^2/(kvec(nn)^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w).*faddeeva(xie2(w),50))+...
                     2*wpi^2/(kvec(nn)^2*vith^2)*(1 + i*sqrt(pi)*xii(w).*faddeeva(xii(w),50)));
end
    z = abs(iadisp(wc,kvec(nn)));
    min(min(z))
    [q r] = find(min(min(z))==z);
    wia(nn) = wrv(r);
    wga(nn) = wiv(q);
    
    wnorm = 1; wpi;
    plot(h(1),kvec(nn)*Ld,wia(nn),'o',kvec(nn)*Ld,wga(nn)/wnorm,'x')
    errorbar(h(1),kvec(nn)*Ld,wrv(fix(nf/2))/wnorm,frange/wnorm)
    errorbar(h(1),kvec(nn)*Ld,wiv(fix(nf/2))/wnorm,frange/wnorm)
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

pos=find(max(wga)==wga)
kmax = kvec(pos(1))
vph = wia(pos(1))/kvec(pos(1))
vphnorm = vph(1)/veth1
maxfreq = wia(pos(1))/wpi
maxgamma = wga(pos(1))/wpi

if 0
%%

plot(kvec*Ld,wia/wpi,kvec*Ld,wga/wpi)
plot(kvec*Ld,wia/wpe,kvec*Ld,wga/wpe)
xlabel('k\lambda_D')
ylabel('\omega_{max}/\omega_{pi}, \gamma_{max}/\omega_{pi}')
distrStr = {['R=' num2str(ne2/ni,'%.2f')],...
            ['Ti=' num2str(Ti,'%.0f')],...
            ['Te1=' num2str(Te1,'%.0f')],...
            ['Te2=' num2str(Te2,'%.0f')],...
            ['S=' num2str(vd2/veth1,'%.2f')],...
            }; 
text(min(get(gca,'xlim')),max(get(gca,'ylim')),distrStr,'horizontalalignment','left','verticalalignment','top')

%find(min(abs(iadisp(wc))))
end