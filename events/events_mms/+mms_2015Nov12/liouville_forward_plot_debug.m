% liouville_forward_plot_debug
% find max f, for scaling color to f.
fmax = 0;
all_f = [];
for iParticle= 1:nParticles 
  all_f = [all_f; saveParticle{iParticle}.f];
  if saveParticle{iParticle}.f >fmax,
    fmax = saveParticle{iParticle}.f;
  end
end
fmax = 0.5e-25;

%%
hca = axes;
holdnow = 1;
for iParticle= 1:10:nParticles 
  if holdnow
    hold(hca,'on')
    holdnow = 0;
    hca.XLabel.String = 'L (km)';
    hca.YLabel.String = 'M (km)';
    hca.ZLabel.String = 'N (km)';
  end
  thisParticle = saveParticle{iParticle};
  thisX = thisParticle.r(:,1)*1e-3; % L
  thisY = thisParticle.r(:,2)*1e-3; % M
  thisZ = thisParticle.r(:,3)*1e-3; % N
  if 0 % colormap according to f
    ncmap = 64;
    cmap = parula(64); % parula
    ncmap = size(cmap,1); % 64   
    colorind = ceil(ncmap*thisParticle.f/fmax);
    if colorind == 0
      thisParticle.f
      plotcolor = [1 1 1];
    elseif colorind > ncmap
      colorind = 64;
      plotcolor = cmap(fix(colorind),:);
    else
      plotcolor = cmap(fix(colorind),:);
    end
  else % colormap according to ending position
    if thisZ(end)<0;
      plotcolor = [0.1000    0.4000    1.0000];
    else
      plotcolor = [0.9500    0.7000         0];
    end
  end
  plot3(hca,thisX,thisY,thisZ,'color',plotcolor)
  pause(0.1)
end
hold(hca,'off')