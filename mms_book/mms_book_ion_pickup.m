% Ion pick-up
vExB = 1;


f = @(vx,vy,vz,ux,uy,uz,vtx,vty,vtz) exp(-((vx-ux)/vtx).^2 - ((vy-uy)/vty).^2 - ((vz-uz)/vtz).^2);


vmax = 2.2;
vx = linspace(-vmax,vmax,151);
vy = linspace(-vmax,vmax,152);
vz = linspace(-vmax,vmax,153);

[VX,VY,VZ] = ndgrid(vx,vy,vz);


hp = findobj(gcf,'type','patch'); delete(hp)


colors = pic_colors('matlab');

h = setup_subplots(1,1);
isub = 1;


hca = h(isub); isub = isub + 1;
if 1 % constant energy surface
  vtx = 1; vty = 1; vtz = 1;
  ux = vExB; uy = 0; uz = 0;
  F = f(VX,VY,VZ,ux,uy,uz,vtx,vty,vtz);
  Flev = f(1,0,0,0,0,0,vtx,vty,vtz);
  s = isosurface(VX,VY,VZ,F,Flev);
  p = patch(hca,'Faces',s.faces,'Vertices',s.vertices,'FaceAlpha',0.99);
  p.EdgeColor = 'none';
  p.FaceColor = [0.2 0.8 0.3];
  p.FaceColor = [0.2 0.2 0.8];
  p.FaceColor = [0 0 0]+0.5;
end

drawnow
if 1 % cold accelerated population
  %hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  ux = [0 vExB-sqrt(1/3) vExB+sqrt(1/3) 2];
  uy = [0 sqrt(1/3) sqrt(1/3) 0];
  uz = [0 -sqrt(1/3) -sqrt(1/3) 0];
  vtx = 0.1; vty = 0.1; vtz = 0.1;
  Flev = f(vtx,vty,vtz,0,0,0,vtx,vty,vtz);
  for ipop = 1:numel(ux)
    F = f(VX,VY,VZ,ux(ipop),uy(ipop),uz(ipop),vtx,vty,vtz);  
    s = isosurface(VX,VY,VZ,F,Flev);
    p = patch(hca,'Faces',s.faces,'Vertices',s.vertices,'FaceAlpha',0.7);
    p.EdgeColor = 'none';
    %p.FaceColor = [0.7 0.3 0.2];
    p.FaceColor = colors(ipop,:);
  end
  hold(hca,'off')
end


hca.XLim = vmax*[-1 1];
hca.YLim = vmax*[-1 1];
hca.ZLim = vmax*[-1 1];
axis(hca,'square')

lighting(hca,'gouraud')
hlight = findobj(hca,'Type','Light');
if isempty(hlight)
  hlight = camlight(hca);
end


hca.XLabel.String = 'v_{L}';
hca.YLabel.String = 'v_{M}';
hca.ZLabel.String = 'v_{N}';
%%
        %s = isocaps(VX,VY,VZ,F,Flev);
        %p = patch(ax,Faces=s.faces,Vertices=s.vertices);
        p = patch(ax,'Faces',s.faces,'Vertices',s.vertices);
        hps(isurf) = p;

        % Default formatting
        p.EdgeColor = 'none';
        p.FaceColor = colors(mod(isurf-1,size(colors,1))+1,:).^0.5; % cycle colors
        p.FaceAlpha = faceAlpha(isurf);
      
      %hold(ax,'off') % not with patch objects

      % might be leftover from previous plottings, only add if no lights
      hlight = findobj(ax,'Type','Light');
      if isempty(hlight)
        hlight = camlight(ax);
      end

      lighting(ax,'gouraud')

      view(ax,[2 2 1])
      %view(ax,[0 0 1])

      all_handles.Patch = hps;
      all_handles.Light = hlight;

      axis(ax,'square')
      axis(ax,'equal')