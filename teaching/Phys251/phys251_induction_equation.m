

L = 1;
B0 = 1;
dt = 1;
B = @(z,L) B0*tanh(z/L);

zvec = L*linspace(-3,3,100);
B1 = B(zvec,L);



fontsize = 12;

nrows = 2;
ncols = 1;
h = gobjects([nrows,ncols]);
ipanel = 1;
for irow = 1:nrows
  for icol = 1:ncols
    h(irow,icol) = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
  end 
end
isub = 1;


if 1 % diffusion
  etamu = 0.2;
  v = 0;
  Bdiff = etamu*gradient(gradient(B1,zvec),zvec);
  Bconv = gradient(B1*v,zvec);
  B2 = B1 + Bdiff + Bconv;

  hca = h(isub); isub = isub + 1;
  plot(hca,zvec,B1,zvec,Bdiff,zvec,Bconv,zvec,B2)
  %legend(hca,{'B_0','(\eta/mu_0)\nabla^2 B','\nabla\times(v\times B)','B = B_0 + \Delta{t} \partial_t{B}'},...
  %  'box','off','location','best')
  irf_legend(hca,{'B_0','(\eta/mu_0)\nabla^2 B','\nabla\times(v\times B)','B = B_0 + \Delta{t} \partial_t{B}'}',...
    [0.98 0.13],'fontweight','bold','fontsize',12)
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  hca.Title.String = 'No convection, v = 0';
end

if 1 % convection
  etamu = 0.0;
  v = 0.1;
  dL = v*dt;
  Bdiff = etamu*gradient(gradient(B1,zvec),zvec);
  Bconv = gradient(B1*v,zvec);
  B2 = B1 + Bdiff + Bconv;
  hca = h(isub); isub = isub + 1;
  plot(hca,zvec,B1,zvec,Bdiff,zvec,Bconv,zvec+0*dL,B2)
  
  %legend(hca,{'B_0','(\eta/mu_0)\nabla^2 B','\nabla\times(v\times B)','B = B_0 + \Delta{t} \partial_t{B}'},...
  %  'box','off','location','best')
  irf_legend(hca,{'B_0','(\eta/mu_0)\nabla^2 B','\nabla\times(v\times B)','B = B_0 + \Delta{t} \partial_t{B}'}',...
    [0.98 0.13],'fontweight','bold','fontsize',12)
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  hca.Title.String = 'No diffusion, \eta = 0';
end

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 12;',1:numel(h))