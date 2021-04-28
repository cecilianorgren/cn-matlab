phi0 = -1;
tau = 1;
xau = 1;
phi = @(x,t) (t/tau).*phi0.*exp(-(x./xau).^2);
ex = @(x,t) 2.*x/(tau^2).*(t/tau).*phi0.*exp(-(x./xau).^2);

nx = 100; xmax = 3;
nt = 100;
t = linspace(0,tau,nt);
x = linspace(-xmax,xmax,nx);
xt = linspace(-xmax,0,nt); 

[X,T] = meshgrid(x,t);
[XT,Txt] = meshgrid(xt,t);
PHI = phi(X,T);
EX = ex(X,T);

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,t,x,PHI')
  hca.XLabel.String = 't';
  hca.YLabel.String = 'x';
  hca.YDir = 'normal';
  hold(hca,'on')
  plot(hca,t,xt,'k')
  hold(hca,'off')
end 
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,t,x,EX')
  hca.XLabel.String = 't';
  hca.YLabel.String = 'x';
  hca.YDir = 'normal';
  hold(hca,'on')
  plot(hca,t,xt,'k')
  hold(hca,'off')
end 
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,t,(phi(xt,t)),t,(ex(xt,t)));
  hca.XLabel.String = 't';
  hca.YLabel.String = 'x';
  legend(hca,{'\phi','E_x=-\partial_x\phi'},'location','best')

end
