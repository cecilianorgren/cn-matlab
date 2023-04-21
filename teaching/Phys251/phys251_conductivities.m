mime = 1836;
n = 1;

cH = @(n,vin,ven,wci,wce) n.*(vin./wci/((vin./wci)^2+1) + ven./wce/((ven./wce)^2+1));

cPs = @(n,xs) n.*(xs./(xs.^2+1));
cP = @(n,xi,xe) cPs(n,xi) + cPs(n,xe);
cHs = @(n,xs) n.*(1./(xs.^2+1));
cH = @(n,xi,xe) cHs(n,xi) - cHs(n,xe);

vin = logspace(-2,log10(4),200); % ~vi
%vin = logspace(log10(1/4),log10(0.001),200); % ~vi
%vin = logspace(-10,log10(4),200);
vevi = 1;
ven = vevi*vin; % ~ve
wci = 1;
wce = mime;
xi = vin/wci;
xe = ven/wce;

%xi = 1./xi;
%xe = 1./xe;

n = xi;
n = xi*0 + 1;



fontsize = 12;

nrows = 4;
ncols = 1;
h = gobjects([nrows,ncols]);
ipanel = 1;
for irow = 1:nrows
  for icol = 1:ncols
    h(irow,icol) = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
  end 
end
isub = 1;


if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,xi,cP(n,xi,xe))
end
if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,xi,cH(n,xi,xe))
end

if 1 % density profile
  hca = h(isub); isub = isub + 1;
  plot(hca,xi,n)
  %legend(hca,{'\sigma_P','\sigma_H'},'location','best','box','off')
  hca.XLabel.String = '\nu_{ni}/\Omega_i';
  hca.YLabel.String = 'n_n';
  irf_legend(hca,{['\nu_{ne}= ' num2str(vevi) ' \nu_{ni}']},[0.3 0.98],'color','k')
end
if 1 % cP, cH
  hca = h(isub); isub = isub + 1;
  plot(hca,xi,cP(n,xi,xe),xi,cH(n,xi,xe))
  legend(hca,{'\sigma_P','\sigma_H'},'location','best','box','off')
  hca.XLabel.String = '\nu_{ni}/\Omega_i';
  irf_legend(hca,{['\nu_{ne}= ' num2str(vevi) ' \nu_{ni}']},[0.3 0.98],'color','k')
end
if 1 % cPi, cPe
  hca = h(isub); isub = isub + 1;  
  plot(hca,xi,cPs(n,xi),xi,cPs(n,xe))
  legend(hca,{'\sigma_{Pi}‚','\sigma_{Pe}'},'location','best','box','off')
  hca.XLabel.String = '\nu_{ni}/\Omega_i';
  irf_legend(hca,{['\nu_{ne}= ' num2str(vevi) ' \nu_{ni}']},[0.3 0.98],'color','k')
end

if 1 % cHi, cHe
  hca = h(isub); isub = isub + 1;  
  plot(hca,xi,cHs(n,xi),xi,cHs(n,xe))
  legend(hca,{'\sigma_{Hi}‚','\sigma_{He}'},'location','best','box','off')
  hca.XLabel.String = '\nu_{ni}/\Omega_i';  
  irf_legend(hca,{['\nu_{ne}= ' num2str(vevi) ' \nu_{ni}']},[0.3 0.98],'color','k')
end



if 0
  hca = h(isub); isub = isub + 1;
  %vevi = vevi/40;
  plot(hca,xi,cP(n,xi,xe),xi,cH(n,xi,xe))
  legend(hca,{'\sigma_P','\sigma_H'},'location','best','box','off')
  hca.XLabel.String = '\nu_{ni}/\Omega_i';
  irf_legend(hca,{['\nu_{ne}= ' num2str(vevi) ' \nu_{ni}']},[0.3 0.98],'color','k')
end












