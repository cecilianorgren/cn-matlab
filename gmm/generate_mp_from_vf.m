function x = generate_mp_from_vf(xe,f,p,Ntot)

%rng(1)
%xc = xe(1:end-1) + 0.5*diff(xe);
nBin = numel(xe)-1;
%fgrid = f(xc);
NBin = round(Ntot*f/sum(f));

x = cell(1,nBin);
for iBin = 1:nBin
  r = rand(NBin(iBin),1);
  x1 = xe(iBin);
  x2 = xe(iBin+1);
  x{iBin} = (x1.^p + r*(x2.^p-x1.^p)).^(1./p);
end