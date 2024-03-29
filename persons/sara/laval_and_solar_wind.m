dvdr = @(v,dAdr,vs) v*(-dAdr)./(1-v.^2/vs^2);

vs = 1;
v = linspace(0,3,1000);
dAdr = -1; % converging

plot(v/vs,dvdr(v,dAdr,vs),'.',v/vs,dvdr(v,-dAdr,vs),'.',v/vs,dvdr(v,0,vs),'.')
legend('converging','diverging')
grid on

%%

fun_r = @(x,x0,rmin,rmax) rmax - rmin*exp(-(x-x0).^2);
fun_A = @(x,x0,rmin,rmax) pi*fun_r(x,x0,rmin,rmax).^2;

xvec = linspace(0,2);
dx = xvec(2) - xvec(1);
x0 = 1;
rmin = 1;
rmax = 2;


nrows = 3;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,fun_r(xvec,x0,rmin,rmax))
hca.YLim = [0 rmax];

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,fun_A(xvec,x0,rmin,rmax))
hca.YLim(1) = 0;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec(2:end),diff(fun_A(xvec,x0,rmin,rmax))/dx)


%% Syms, Venturi/Venturi-Laval
% Define variables
syms x v(x) t
Dv = diff(v,x);

% Define constants
rmax = 2;
rmin = 1.8506; % v0 = 0.2vs
rmin = 0.5;

vs = 1;
x0 = 3;
xvec = linspace(0,2*x0,1000);

% Specify box
r = rmax - (rmax-rmin)*exp(-(x-x0)^2);
r = 0.1*x/(rmax-rmin);
A = pi*r^2;
dAdx = gradient(A,x);

% Specify and solve ode
ode = Dv == v*(-dAdx/A)/(1-v^2/vs^2);
cond1 = v(x0) == 1.0*vs;
cond2 = Dv(x0) == 0;
conds = [cond1];
%cond = r(0) == 0;
vSolx = dsolve(ode,conds);
dvdx = gradient(vSolx,x);

% Numerical
%rhs = matlabFunction(v*(-dAdx/A)/(1-v^2/vs^2));
%v_num = ode45(rhs,[0 x0],[0.1 0.1]);


% Make functions
mf_r = matlabFunction(r);
mf_A = matlabFunction(A);
mf_dAdx = matlabFunction(dAdx);
mf_v = matlabFunction(vSolx);
mf_dvdx = matlabFunction(dvdx);

% Plot
nrows = 5;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_r(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'r';
hca.YLim = [0 rmax];

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_A(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'A';
hca.YLim(1) = 0;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_dAdx(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'dA/dx';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,real(mf_v(xvec)),xvec,imag(mf_v(xvec)))
hca.XLabel.String = 'x';
hca.YLabel.String = 'v';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
%plot(hca,mf_v(xvec),real(mf_v(xvec)),xvec,imag(mf_v(xvec)))
plot(hca,mf_v(xvec)/vs,mf_dvdx(xvec))
hca.XLabel.String = 'v/v_s';
hca.YLabel.String = 'dv/dx';

%% Numerical, Venturi/Venturi-Laval
% Define variables
syms x v(x)
Dv = diff(v,x);

% Define constants
rmax = 2;
rmin = 1.8506; % v0 = 0.2vs
rmin = 1.5;

vs = 1;
x0 = 3;
xvec = linspace(0,2*x0,1000);

% Specify box
r = rmax - (rmax-rmin)*exp(-(x-x0)^2);
%r = 0.1*x/(rmax-rmin);
A = pi*r^2;
dAdx = gradient(A,x);

% Specify and solve ode
ode = Dv == v*(-dAdx/A)/(1-v^2/vs^2);
cond1 = v(x0) == 1*vs;
cond2 = Dv(x0) == 0;
conds = [cond1];
%cond = r(0) == 0;
vSolx = dsolve(ode,conds);
dvdx = gradient(vSolx,x);

% Make functions
mf_r = matlabFunction(r);
mf_A = matlabFunction(A);
mf_dAdx = matlabFunction(dAdx);
mf_v = matlabFunction(vSolx);
mf_dvdx = matlabFunction(dvdx);

% Plot
nrows = 5;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_r(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'r';
hca.YLim = [0 rmax];

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_A(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'A';
hca.YLim(1) = 0;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,mf_dAdx(xvec))
hca.XLabel.String = 'x';
hca.YLabel.String = 'dA/dx';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,xvec,real(mf_v(xvec)),xvec,imag(mf_v(xvec)))
hca.XLabel.String = 'x';
hca.YLabel.String = 'v';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
%plot(hca,mf_v(xvec),real(mf_v(xvec)),xvec,imag(mf_v(xvec)))
plot(hca,mf_v(xvec)/vs,mf_dvdx(xvec))
hca.XLabel.String = 'v/v_s';
hca.YLabel.String = 'dv/dx';

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
