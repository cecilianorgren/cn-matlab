function out = tearing_bx(xvec,kvec,a,C)

bx = xvec*0;

Amat = C./(2*kvec.*a).*exp(-2*kvec.*a);
Bmat = C./(2*kvec.*a).*(2*kvec.*a-1);

for ip = 1:numel(xvec)
  k = kvec(ip);
  x = xvec(ip);
  A = Amat(ip);
  B = Bmat(ip);
  if x > a
    bx(ip) = C*exp(-k.*x);
  elseif x < a && x > 0
    bx(ip) = A*exp(k.*x) + B*exp(-k.*x);
  elseif x < 0 && x > -a
    bx(ip) = A*exp(-k.*x) + B*exp(k.*x);
  elseif x < -a
    bx(ip) = C*exp(k.*x);
  end
end

out = bx;
