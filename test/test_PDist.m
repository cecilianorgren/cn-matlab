pdist = ePDist3.tlim(times(1)+0.015*[-1 1]);
pdist = iPDist3.tlim(times(1)+0.075*[-1 1]);

pdist = ePDist3(1).elim([100 Inf]);

L = [1 1 0]; L = L/norm(L);
N = cross(L,M); N = N/norm(N);
M = cross(N,L); M = M/norm(M);

%L = [1 0 0]; 
%M = [0 1 0]; 
%N = [0 0 1]; 

nMC = 1000;
aa = pdist.shift(10000*M,nMC,[L;M;N],'mms');
aa0 = pdist.shift(000*M,nMC,[L;M;N],'mms');
bb = pdist.convertto('s^3/m^6');

vg = -70000:2000:70000;
%vg = -3000:200:3000;


aar  =  aa.reduce('2D',[1 0 0],[0 1 0],'vg',vg,'nMC',100);
aa0r = aa0.reduce('2D',[1 0 0],[0 1 0],'vg',vg,'nMC',100);
bbr  =  bb.reduce('2D',L,M,'vg',vg,'nMC',100);

% Plot, reduced dists
h = gobjects(0);
nrows = 1;
ncols = 3;
ip = 0;
for irows = 1:nrows
  for icols = 1:ncols
    ip = ip + 1;
    h(ip) = subplot(nrows,ncols,ip);
  end
end
h = reshape(h,icols,irows);
h = transpose(h);
h = h(:);

isub = 1;

if 0 % non-formatted
  hca = h(isub); isub = isub + 1;
  imagesc(hca,log10(squeeze(aar.data)))
  hcb = colorbar(hca);
  hca.Title.String = 'non-zero shift';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,log10(squeeze(aa0r.data)))
  hcb = colorbar(hca);
  hca.Title.String = 'zero shift';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,log10(squeeze(bbr.data)))
  hcb = colorbar(hca);
  hca.Title.String = 'original';
end
if 1 % formatted
  hca = h(isub); isub = isub + 1;
  aar.plot_plane(hca)
  hca.Title.String = 'non-zero shift';

  hca = h(isub); isub = isub + 1;
  aa0r.plot_plane(hca)
  hca.Title.String = 'zero shift';

  hca = h(isub); isub = isub + 1;
  bbr.plot_plane(hca)
  hca.Title.String = 'original';
end
hlinks2 = linkprop(h,{'CLim'});
c_eval('axis(h(?),''square'')',1:numel(h))
c_eval('axis(h(?),''equal'')',1:numel(h))
c_eval('h(?).XLim = vg([1 end])*1e-3; h(?).YLim = vg([1 end])*1e-3;',1:numel(h))

%% Plot, skymaps
pdist = iPDist3(1);
nMC = 1000;
aa0 = pdist.shift(000*[0 1 0],nMC,[1 0 0;0 1 0; 0 0 1],'mms');
bb = pdist.convertto('s^3/m^6');

h = gobjects(0);
nrows = 2;
ncols = 1;
ip = 0;
for irows = 1:nrows
  for icols = 1:ncols
    ip = ip + 1;
    h(ip) = subplot(nrows,ncols,ip);
  end
end

iEnergy = 22;

isub = 1;

hca = h(isub); isub = isub + 1;
imagesc(hca,squeeze(aa0.data(1,iEnergy,:,:)))
hcb = colorbar(hca);

hca = h(isub); isub = isub + 1;
imagesc(hca,squeeze(bb.data(1,iEnergy,:,:)))
hcb = colorbar(hca);

hlinks2 = linkprop(h,{'CLim'});

%%
pdist = iPDist3(1);
nMC = 1000;
aa0 = pdist.shift(000*[0 1 0],nMC,[1 0 0;0 1 0; 0 0 1],'mms');
aa0 = aa0;
bb = pdist.convertto('s^3/m^6');
nrows = 2;
ncols = 1;
ip = 0;
for irows = 1:nrows
 for icols = 1:ncols
  ip = ip + 1;
  h(ip) = subplot(nrows,ncols,ip);
 end
end
iEnergy = 22;
isub = 1;
hca = h(isub); isub = isub + 1;
imagesc(hca,squeeze(aa0.data(1,iEnergy,:,:)))
hcb = colorbar(hca);
hca = h(isub); isub = isub + 1;
imagesc(hca,squeeze(bb.data(1,iEnergy,:,:)))
hcb = colorbar(hca);
caxis(h(2),caxis(h(1)))