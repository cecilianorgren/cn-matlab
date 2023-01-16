% illustrate_magnetic_reconnection

doVideo = 1;
doGif = 1;
fileName = 'illustration_magnetic_reconnection';

a = 5;
b = 1;
x = a*linspace(-10,10,200);
y = b*linspace(-10,10,100);
z = linspace(-10,10,5);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);
x_xline = x;
y_xline = x*b/a;

Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY0 = Ay(X,Y);

%[FX,FY] = gradient(AY,dx,dy);
%Bx = -FX;
%By = FY;

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

hca = subplot(1,1,1);
t = 0:30;
Astep = 20;
dA = Astep/numel(t);
AYlev0 = -100:Astep:(100 + Astep);

% Initiate
if doVideo
  vidObj = VideoWriter([printpath fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end
     

for it = 1:numel(t)
  % Draw separatrix
  plot(hca,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot(hca,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])

  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot(hca,sx(sy>=0),sy(sy>=0),'color',[1 0 0],'linewidth',2)
    plot(hca,sx(sy<=0),sy(sy<=0),'color',[0 0 1],'linewidth',2)
   

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
