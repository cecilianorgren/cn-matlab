no02m = PIC('/Users/cecilianorgren/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');


%% Move of ion temperature

pic = no02m.twpelim([15000 24000]).xlim([40 160]).zlim([-10 10]);

%pic.movie({'ti'},'A',1,'clim',{[0 0.6*0.99]},'cmap',{pic_colors('thermal')},'filename',[printpath 'no02m_Ti'])
%pic.movie({'ti'},'A',1,'clim',{[0 0.6*0.99]},'cmap',{pic_colors('thermal')},'filename',[printpath 'no02m_Ti_2'])

pic.movie({'ti','te'}','A',1,'clim',{[0 0.6*0.99],[0 0.6*0.99]},'cmap',{pic_colors('thermal'),pic_colors('thermal')},'filename',[printpath 'no02m_TiTe_2'])


%%
varstrs = {'ti','te','Babs'};
cmaps = {pic_colors('thermal'),pic_colors('thermal'),pic_colors('thermal')};

no02m.twpelim(20000).xlim([60 140]).zlim([-5 5]).plot_map(varstrs','A',0.2,'cmaps',cmaps)

%% Ti, Te histogram
varstrs = {'ti','te','Babs'};
cmaps = {pic_colors('thermal'),pic_colors('thermal'),pic_colors('thermal')};

pic = no02m.twpelim(22000).xlim([60 140]).zlim([-5 5]);
ti = pic.t([3 5]);
te = pic.t([4 6]);
%ti = pic.t([1]);
%te = pic.t([2]);
Ay = pic.A;

[Ainds,Avals] = saddle(Ay,'sort');
Ay_xline = Avals(1);
Ay_presheet = 10;

for iA = 0
%iA = -1;
Amin = Ay_xline + 0.5*(iA-1);
Amax = Ay_xline + 0.5*iA;

ti(Ay<Amin) = NaN;
te(Ay<Amin) = NaN;

ti(Ay>Amax) = NaN;
te(Ay>Amax) = NaN;

ti(ti<0) = 0;
te(te<0) = 0;

Alevels = Ay_xline + -10:0.5:10;


ti_edges = linspace(0,0.7,100);
te_edges = linspace(0,0.4,100);



%ti_edges = linspace(0,0.15,100); % for pressure
%te_edges = linspace(0,0.15,100);

%histcn_plot([ti(:), te(:)], t_edges,t_edges)
[count edges mid loc] = histcn([ti(:), te(:)], ti_edges,te_edges(2:end));

hca = subplot(2,4,[1 2]);
pcolor(hca,pic.xi,pic.zi,ti')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = 'Ion temperature';
hold(hca,'on')
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')

hca = subplot(2,4,[5 6]);
pcolor(hca,pic.xi,pic.zi,te')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = 'Electron temperature';
hold(hca,'on')
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')


hca = subplot(2,4,[3 4 7 8]);
hca.Position = [0.57 0.15 0.3 0.7];

pcolor(hca,mid{1},mid{2},log10(count)')
shading(hca,'flat')
%hca.CLim = prctile(count(:),[0 99.7]);
hcb = colorbar(hca);
hcb.YLabel.String = 'log_{10} Counts';
hca.Box = 'on';
hca.Layer = 'top';
hca.FontSize = 16;
hca.XLabel.String = 'Ion temperature';
hca.YLabel.String = 'Electron temperature';

colormap(pic_colors('thermal'))

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 16;',1:numel(h))
drawnow
pause(0.1)
cn.print(sprintf('no02m_twpe_22000_ne_te_exhaust_A%g',iA))
end
%% general 2D histogram

pic = no02m.twpelim(22000).xlim([60 140]).zlim([-5 5]);
varstr1 = {'ne'};
varstr2 = {'te'};
var1 = pic.n([4 6]); 
var2 = pic.t([4 6]); 

var1(var1<0) = 0;
var2(var2<0) = 0;

var1_min = min(var1(:));
var1_max = max(var1(:));
var2_min = min(var2(:));
var2_max = max(var2(:));
var1_min = min(var1(:));
var1_max = prctile(var1(:),99);
var2_min = min(var2(:));
var2_max = prctile(var2(:),99);
var1_edges = linspace(var1_min,var1_max,100);
var1_edges = linspace(var2_min,var2_max,100);

Ay = pic.A;

[Ainds,Avals] = saddle(Ay,'sort');
Ay_xline = Avals(1);
Ay_presheet = 10;

for iA = 10
%iA = -1;
Amin = Ay_xline + 0.5*(iA-1);
Amax = Ay_xline + 0.5*iA;

var1(Ay<Amin) = NaN;
var2(Ay<Amin) = NaN;

var1(Ay>Amax) = NaN;
var2(Ay>Amax) = NaN;

%var1(ti<0) = 0;
%te(te<0) = 0;

Alevels = Ay_xline + -10:0.5:10;


% p*n^(-gamma) = const
% T = const: gamma = 1 -> p*n^(-1) = p/n = nT/n = T = const.
% adiabatic isotropic: gamma = 5/3 -> p*n^(-5/3) = Tn^(3/3)*n^(-5/3) = T*n^(-2/3) = const.
EqSt = @(n,T,gamma) T.*n.*n.^(-gamma);
gamma = 5/3;

EqSt = @(n,T,gamma) T.*n;
gamma = 5/3;


%histcn_plot([ti(:), te(:)], t_edges,t_edges)
[count edges mid loc] = histcn([var1(:), var2(:)], var1_edges,var1_edges(2:end));

hca = subplot(2,4,[1 2]);
pcolor(hca,pic.xi,pic.zi,var1')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = varstr1;
hold(hca,'on')
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')

hca = subplot(2,4,[5 6]);
pcolor(hca,pic.xi,pic.zi,var2')
shading(hca,'flat')
hca.XLabel.String = 'x (c/\omega_{pi})';
hca.YLabel.String = 'z (c/\omega_{pi})';
hcb = colorbar(hca);
hcb.YLabel.String = varstr1;
hold(hca,'on')
clim = hca.CLim;
contour(hca,pic.xi,pic.zi,Ay',Alevels,'k')
hca.CLim = clim;
hold(hca,'off')


hca = subplot(2,4,[3 4 7 8]);
hca.Position = [0.57 0.15 0.3 0.7];

pcolor(hca,mid{1},mid{2},log10(count)')
shading(hca,'flat')
%hca.CLim = prctile(count(:),[0 99.7]);
hcb = colorbar(hca);
hcb.YLabel.String = 'log_{10} Counts';
hca.Box = 'on';
hca.Layer = 'top';
hca.FontSize = 16;
hca.XLabel.String = varstr1;
hca.YLabel.String = varstr2;
hold(hca,'on')
clim = hca.CLim;
[X,Z] = ndgrid(mid{1},mid{2});
EQST = EqSt(X,Z,gamma);
contour(hca,X,Z,EQST,50,'w')
hca.CLim = clim;
hold(hca,'off')

colormap(pic_colors('thermal'))

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 16;',1:numel(h))
drawnow
pause(0.1)
cn.print(sprintf('no02m_twpe_22000_ne_te_exhaust_A%g_p=const',iA))
end

%% Earth as a plasma laboratory

%% Movie of reconnection, sketch
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
%By = FY;1


colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

hca = subplot(1,1,1);
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
  plot(hca,x_xline,y_xline,'linewidth',4,'linestyle','--','color',color_separatrix)
  hold(hca,'on')
  plot(hca,x_xline,-y_xline,'linewidth',4,'linestyle','--','color',color_separatrix)

  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
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
  
  % Draw magnetic energy density
  [BZ,BX] = gradient(AY,dy,-dx);
  %BZ = gradient(AY,dx);
  UB = (BX.^2 + BZ.^2)/2;
  [X,Y] = ndgrid(x,y);
  pcolor(hca,X,Y,(UB)')
  shading(hca,'interp')
  colormap(pic_colors('blue_red'))
  %colormap(hca,'gray')
  hca.CLim = [0 80];

  hca.Children = circshift(hca.Children,-1);
  

  pause(0.01)
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

%%
N = 200;
[X,Y,Z] = sphere(N);

hca = subplot(1,1,1);
colors = pic_colors('matlab');
surf(hca,X,Y,Z)
shading(hca,'flat')
colormap(colors(6,:))
camlight
axis(hca,'equal')
hca.Visible = 'off';  