% illustrate_magnetic_reconnection

doVideo = 1;
doGif = 1;
fileName = 'illustration_magnetic_reconnection';

a = 5;
b = 1;
x = a*linspace(-10,10,200);
z = b*linspace(-10,10,100);

[X,Z] = meshgrid(x,z);
dx = x(2) - x(1);
dz = z(2) - z(1);

x_xline = x;
z_xline = x*b/a;

Ay = @(x,z) (x/a).^2 - (z/b).^2;
AY0 = Ay(X,Z);

%[FX,FY] = gradient(AY,dx,dy);
%Bx = -FX;
%By = FY;1


colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

hca = subplot(1,1,1);
hca.Position(4) = hca.Position(3)*(b/a);

t = 0:60;
Astep = 20;
dA = Astep/numel(t);
AYlev0 = -100:Astep:(100 + Astep);

linecolor_top = colors(1,:);
linecolor_bot = colors(2,:);
color_background = [0 0 0];
color_separatrix = [0 0 0];

% Initiate
if doVideo
  vidObj = VideoWriter([printpath fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end
     
% Advance and collect frames
for it = 1:numel(t)
  % Draw separatrix
  plot(hca,x_xline,z_xline,'linewidth',4,'linestyle','--','color',color_separatrix)
  %axis(hca,'equal')
  hold(hca,'on')
  plot(hca,x_xline,-z_xline,'linewidth',4,'linestyle','--','color',color_separatrix)

  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,z,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.1:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.1:numel(S(is).Y));%S(is).Y;
    plot(hca,sx(sy>=0),sy(sy>=0),'color',linecolor_top,'linewidth',8)
    plot(hca,sx(sy<=0),sy(sy<=0),'color',linecolor_bot,'linewidth',8)
   

%     if not(holdon)
%       hold(hca,'on')
%       holdon = 1;
%     end
  end
  pause(0.1)
  drawnow
  hold(hca,'off')
  hca.Visible = 'off';
  hca.Position = [0 0 1 1];

  % Collect frames
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = numel(t);
      currentBackgroundColor = get(gcf,'color');
      set(gcf,'color',[1 1 1]);
      drawnow      
      tmp_frame = getframe(gcf);
      %cell_movies{imovie}(itime) = tmp_frame;
      if iframe == 1 % initialize animated gif matrix
        [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
        %map(end+1,:) = get(gcf,'color');
        im_tmp(1,1,1,nframes) = 0;                                                
        all_im = im_tmp;             
      else
        all_im(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,map,'nodither');
      end       
    end    
  end
end

% Finalize
if doVideo
  close(vidObj)
end
if doGif
  imwrite(all_im,map,[printpath fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
end
