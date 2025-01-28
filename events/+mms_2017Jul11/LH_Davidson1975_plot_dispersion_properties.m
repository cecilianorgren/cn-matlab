figure(100)

% Set up plot
nrows = 8;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

isub = 1;
if 1 % b
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,b(kvec,vte,wce))
  hca.YLabel.String = 'b(k,vte,wce)';
end
if 1 % bessel0(b)
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,besseli(0,b(kvec,vte,wce)))
  hca.YLabel.String = 'I_0(b)';
end
if 1 % bessel1(b)
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,besseli(1,b(kvec,vte,wce)))
  hca.YLabel.String = 'I_1(b)';
end
if 1 % bessel1(b)/bessel0(b)
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,besseli(1,b(kvec,vte,wce))./besseli(0,b(kvec,vte,wce)))
  hca.YLabel.String = 'I_1(b)/I_0(b)';
end
if 1 % b*(1-bessel1(b)/bessel0(b))
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,b(kvec,vte,wce).*(1-besseli(1,b(kvec,vte,wce))./besseli(0,b(kvec,vte,wce))))
  hca.YLabel.String = 'b(1-I_1(b)/I_0(b))';
end
if 1 % 1-b*(1-bessel1(b)/bessel0(b))
  hca = h(isub); isub = isub + 1;
  plot(hca,kvec*roe,1-b(kvec,vte,wce).*(1-besseli(1,b(kvec,vte,wce))./besseli(0,b(kvec,vte,wce))))
  hca.YLabel.String = '1-b*(1-I_1(b)/I_0(b))';
end