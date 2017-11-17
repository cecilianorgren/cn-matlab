%% Derivative of besseli
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

s = 0; % order
x = linspace(0,10,1000); % argument

isub = 1;
if 1 % besselj(s,x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x,besselj(s,x))
  hca.YLabel.String = sprintf('J_%g(x)',s);
end
if 1 % besselj(s-1,x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x,besselj(s-1,x))
  hca.YLabel.String = sprintf('J_{%g}(x)',s-1);
end
if 1 % besselj'(s,x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x,(s*besselj(s, x))./x - besselj(s + 1, x),x(1:end-1),diff(besselj(s,x))/(x(2)-x(1)))
  hca.YLabel.String = sprintf('J''_{%g}(x)',s);
end
