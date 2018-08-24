% dntrap_all{ih} = dntrap;
% x_all{ih} = x;
% phi_all{ih} = phi;
% fun_fit_all{ih} = fun_net;

h(1) = subplot(1,2,1); set(h(1),'LineStyleOrder',{'-','--'});
h(2) = subplot(1,2,2); set(h(2),'LineStyleOrder',{'-','--'});

neh = numel(fun_fit_all);

for ih = 1:neh  
  hca = h(1);
  hold(hca,'on')
  plot(hca,x_all{ih},phi_all{ih})
  hold(hca,'off')
  
  
  hca = h(2);
  hold(hca,'on')  
  plot(hca,phi_all{ih},dntrap_all{ih},'.')
  plot(hca,phi_all{ih},fun_fit_all{ih}(phi_all{ih}))
  hold(hca,'off')
end