function varargout = D_streaming(w,k,vt,wp,vd)
c = 299792458;
Nfad = 50;
nsp = numel(vt);

D = 1;
D_sep = zeros(1,nsp);
xi_arg = @(w,k,vt,wp,vd) (w-k*vd)./k./vt;
xi = @(w,k,vt,wp,vd) 2*wp^2./(k.^2)/(vt^2).*(1+xi_arg(w,k,vt,wp,vd).*i*sqrt(pi).*faddeeva(xi_arg(w,k,vt,wp,vd),Nfad));
for isp = 1:nsp
  D = D + xi(w,k,vt(isp),wp(isp),vd(isp));
  D_sep(isp) = xi(w,k,vt(isp),wp(isp),vd(isp));
end

varargout{1} = D;
if nargout == 2
  varargout{2} = D_sep;
end