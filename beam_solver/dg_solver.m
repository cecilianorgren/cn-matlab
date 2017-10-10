function [wia,wga,kmax,maxfreq,maxgamma,vph,residual] = dg_solver(vd1,vd2,veth1,veth2,vith,wpe,wpe1,wpe2,wpi,kvec,doposwi)
% function [wia,wga,kmax,maxfreq,maxgamma,vph,residual] = dg_solver(vd1,vd2,veth1,veth2,vith,wpe,wpe1,wpe2,wpi,kvec,doposwi)

doDiagnostics = 0;

% Initial range of complex frequencies at low k. Linspace range may need to
% be varied depending on mode.
wapprox = kvec(1)*vd2;
wrv = linspace(0,wapprox*2,1501);
wiv = linspace(0,wapprox*2,1501);
[wr,wi] = meshgrid(wrv,wiv);
wc = wr+i*wi;
nk = numel(kvec);
dk = kvec(2)-kvec(1);

% Dispersion relation
xie1 = @(omega)((omega-kvec(1)*vd1)/(kvec(1)*veth1));
xie2 = @(omega)((omega-kvec(1)*vd2)/(kvec(1)*veth2));
xii = @(omega)(omega/(kvec(1)*vith));
iadisp = @(w) (1+2*wpe1^2/(kvec(1)^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w).*faddeeva(xie1(w),50))+...
                 2*wpe2^2/(kvec(1)^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w).*faddeeva(xie2(w),50))+...
                 2*wpi^2/(kvec(1)^2*vith^2)*(1 + i*sqrt(pi)*xii(w).*faddeeva(xii(w),50)));
z = abs(iadisp(wc));
residul(1)=min(min(z));
[q r] = find(min(min(z))==z);
wia = nan(1,nk);
wga = nan(1,nk);
wia(1) = wrv(r);
wga(1) = wiv(q);

gradr = wia(1);
gradi = wiv(1);
    
xie1 = @(omega,k)((omega-k*vd1)/(k*veth1));
xie2 = @(omega,k)((omega-k*vd2)/(k*veth2));
xii = @(omega,k)(omega/(k*vith));
iadisp = @(w,k) (1+2*wpe1^2/(k^2*veth1^2)*(1 + i*sqrt(pi)*xie1(w,k).*faddeeva(xie1(w,k),50))+...
                 2*wpe2^2/(k^2*veth2^2)*(1 + i*sqrt(pi)*xie2(w,k).*faddeeva(xie2(w,k),50))+...
                 2*wpi^2/(k^2*vith^2)*(1 + i*sqrt(pi)*xii(w,k).*faddeeva(xii(w,k),50)));

if doDiagnostics
  h(1) = subplot(1,2,1);
  h(2) = subplot(1,2,2);
end
               
frange = 50;
%frange = 1e-1*wpe;
for nn = 2:length(kvec)
    frange = vd2*dk;
    try
    %linear search for zero point around linearly interpolated guess. frange
    %may need to be changed depending on mode. Large frange decreases accuracy
    wrv = linspace(wia(nn-1)+gradr-frange,wia(nn-1)+gradr+frange,601);
    wiv = linspace(wga(nn-1)+gradi-frange,wga(nn-1)+gradi+frange,601);
    wrv(wrv<0)=[];   
    if doposwi
        if nn < 5; wiv(wiv<0)=[]; end; 
    end
        
    [wr,wi] = meshgrid(wrv,wiv);
    wc = wr+i*wi;
    
    z = abs(iadisp(wc,kvec(nn)));
        
    if doDiagnostics
      pcolor(h(1),wr,wi,real(iadisp(wc,kvec(nn)))); shading(h(1),'flat'); colorbar('peer',h(1))
      pcolor(h(2),wr,wi,imag(iadisp(wc,kvec(nn)))); shading(h(2),'flat'); colorbar('peer',h(2))
      pause;
    end
      
    residual(nn) = min(min(z));
    [q r] = find(min(min(z))==z);
    
    wia(nn) = wrv(r);
    wga(nn) = wiv(q);
    
    %To do: Write better root finder.    
    gradr = wia(nn)-wia(nn-1);
    gradi = wga(nn)-wga(nn-1);
    catch
        continue
    end
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

pos=find(max(wga)==wga);
kmax = kvec(pos(1));
vph = wia(pos(1))/kvec(pos(1));
vphnorm = vph(1)/veth1;
maxfreq = wia(pos(1))/wpi;
maxgamma = wga(pos(1))/wpi;
