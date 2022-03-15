

iDist = 2000;
dist =  iPDist1(iDist);
particles = dist.macroparticles('ntot',100,'skipzero',1);
vdf = dist.reduce('2D',[1 0 0],[0 1 0]);

h = setup_subplots(1,2);
isub = 1;
if 1 % macro particles
  hca = h(isub); isub = isub + 1;
  scatter(hca,particles.vx,particles.vy,'.k')  
  axis(hca,'square')
end
if 1 % distribution
  hca = h(isub); isub = isub + 1;
  pcolor(hca,vdf.depend{1},vdf.depend{2},squeeze(vdf.data)')
  shading(hca,'flat')
  colorbar('peer',hca)
  axis(hca,'square')
end

hlinks = linkprop(h,{'XLim','YLim'});