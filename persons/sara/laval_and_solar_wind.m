dvdr = @(v,dAdr,vs) v*(-dAdr)./(1-v.^2/vs^2);

vs = 1;
v = linspace(0,3,1000);
dAdr = -1; % converging

plot(v/vs,dvdr(v,dAdr,vs),'.',v/vs,dvdr(v,-dAdr,vs),'.')
legend('converging','diverging')
grid on

%% Solar wind solution
vs = 1;
rc = 1;

fun = @(r,v,c) (v/vs).^2 - log((v/vs).^2) - 4*r/rc + 4*log(r/rc) - c;
%fun = @(r,v,c) (v/vs).^2 - log((v/vs).^2) - 4*r/c - 4*log(r/rc) - c;

v = vs*linspace(0,5,500);
r = rc*linspace(0,5,500);
v = vs*logspace(-2,log10(3),500);
r = rc*logspace(-2,log10(5),500);
[R,V] = meshgrid(r,v);


cvec = -9:2:3;
nrows = numel(cvec);
ncols = 1;

for ic = 1:numel(cvec)
  c = cvec(ic);
  hca = subplot(nrows,ncols,ic);
  vpow = 1;
  pcolor(hca,R,V.^vpow,fun(R,V,c))
  hca.CLim = [-5 5];
  shading(hca,'flat')
  hca.XLabel.String = 'r/r_c';
  hca.YLabel.String = sprintf('(v/v_s)^%.0f',vpow);
  colormap(pic_colors('blue_red'))
  hca.Title.String = sprintf('C = %4.0f',c);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = {'f(r,v)'};
end
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
h(1).Title.String = {'f(r,v) = (v/v_s)^2 - log((v/v_s)^2)-4r/r_c + 4log(r/r_c) - C', h(1).Title.String};

%% contourf plot isntead
cvec = -9:2:3;
nrows = 1;
ncols = 1;
legs = {};
colors = pic_colors('matlab');
hca = subplot(nrows,ncols,1);

for ic = 1:numel(cvec)
  c = cvec(ic);    
  vpow = 1;  
  hc = contourf(hca,R,V.^vpow,fun(R,V,c),[0 0],'color',colors(ic,:),'linewidth',1,'tag',num2str(c));  
  shading(hca,'flat')
  hca.XLabel.String = 'r/r_c';
  hca.YLabel.String = sprintf('(v/v_s)^%.0f',vpow);
  %colormap(pic_colors('blue_red'))
  legs{ic} = sprintf('C = %4.0f',c);
  %hcb = colorbar('peer',hca);
  %hcb.YLabel.String = {'f(r,v)'};
  if ic == 1, hold(hca,'on'); end
  if ic == numel(cvec), hold(hca,'off'); end
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
h(1).Title.String = {'f(r,v) = (v/v_s)^2 - log((v/v_s)^2)-4r/r_c + 4log(r/r_c) - C = 0', h(1).Title.String};
legend(hca,legs,'location','east')
grid(hca,'on')
hca.Position(2) = 0.15;
hca.Position(4) = 0.70;


%% fminsearch
c = -3;
v = 2;
FUN = @(r)fun(r,v,c);
X = fminsearch(FUN,1);
