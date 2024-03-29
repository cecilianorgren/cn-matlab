% Data
Lz = 1;
Lx = [0.5,1,2,5,6,7,8,10,15,20,40];
R = [0.22602171,0.24566273,0.25104229,0.20345172,0.18692082,0.17226913,0.15952006,0.14030094,0.09496913,0.05082819,0.03705537];
% Fit function
xx = linspace(0,40,1000);
Lf = 0.18;
lexp = 0.4;
f0 = 1.8;
f = @(Lx) f0*Lx.*exp(-(Lx/Lf).^lexp);
% Lf = 3;
% lexp = 0.4;
% f0 = 1.5;
% f = @(Lx) f0*(Lx).^(lexp-1)./Lf*lexp.*exp(-(Lx/Lf).^lexp);
funstr = '';
model = 7;
switch model
  case 0
    Lf = 0.17; % b
    lexp = 0.4; % c
    f0 = 1.7; % a
    exp0 = 0.02; % zero level
    f = @(Lx) f0*Lx.*exp(-((Lx/Lf)).^lexp)+exp0;    
  case 1
    Lf = 0.18;
    lexp = 0.4;
    f0 = 1.6;
    f = @(Lx) f0*Lx.*exp(-(Lx/Lf).^lexp);
  case 2 % Judit
    b = 0.026; % b
    a = 4.838; % a
    c = 0.299; % c
    f = @(Lx) a*Lx.*exp(-(Lx/b).^c);
  case 3 
    b = 1.6; % b
    a = 0.4; % a
    c = 0.55; % c
    E0 = 0.03;
    f = @(Lx) a*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^c) + E0;
  case 4 
    b = 8; % b
    a = 0.3; % a
    c = 0.2; % c
    E0 = 0.03;
    f = @(Lx) a*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^(c+1)) + E0;
    funstr = {'f(L_x) = a(c+1)(L_x/b)^c exp[-(L_x/b)^{c+1}]+f_{0,inf}';...
      sprintf('a = %g, b = %g, c = %g, f_{0,inf} = %g',a,b,c,E0)};
  case 5 
    b = 8; % b
    a = 2.4; % a
    c = 0.2; % c
    E0 = 0.03;
    f = @(Lx) a/b*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^(c+1)) + E0;
    funstr = {'f(L_x) = (a/b)(c+1)(L_x/b)^c exp[-(L_x/b)^{c+1}]+f_{0,inf}';...
      sprintf('a = %g, b = %g, c = %g, f_{0,inf} = %g',a,b,c,E0)};
  case 6 
    b = 8; % b
    a = 2.4*8; % a
    c = 0.18; % c
    E0 = 0.03;
    f = @(Lx) a/b/b*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^(c+1)) + E0;
    funstr = {'f(L_x) = (a/b)(c+1)(L_x/b)^c exp[-(L_x/b)^{c+1}]+f_{0,inf}';...
      sprintf('a = %g, b = %g, c = %g, f_{0,inf} = %g',a,b,c,E0)};
  case 7 
    b = 10; % b
    a = 2.4*8; % a
    c = 0.1; % c
    E0 = 0.03;
    e = 0.4;
    f = @(Lx) a/b/b*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^(c+1)+e) + E0;
    funstr = {'f(L_x) = (a/b)(c+1)(L_x/b)^c exp[-(L_x/b)^{c+1} + e]+f_{0,inf}';...
      sprintf('a = %g, b = %g, c = %g, f_{0,inf} = %g, e = %g',a,b,c,E0,e)};
end
    

%f = @(Lx) f0*Lx.*(exp(-((Lx/Lf)).^lexp)+exp0./Lx);


% Judits fit

%f = @(Lx) f0*Lx.*(exp(-((Lx/Lf)).^lexp)+exp0./Lx);
L = sqrt(Lx.^2 + Lz.^2);

Lz = 1;
%f = @(Lx) Lz./(Lx+3).^0.5-0.15;
%f = @(Lx) (Lz./(Lx+5)).^0.5-0.12;
%f = @(Lx) (Lz./(Lx)).^0.5-0.15;
%f = @(Lx) atan((Lz)./(Lx)+0.02).^0.95;
%f = @(Lx) atan(2*(Lz)./(Lx)+1.00).^1.0;

%Lz/Lx

plot(Lx,R,'*',xx,f(xx),'-')
%plot(L,R,'*')
%semilogx(Lx,R,'*',xx,f(xx),'-')
hca = gca;
hca.YLim = [0 0.3];
grid on;

irf_legend(hca,funstr,[0.98 0.98],'color',[0 0 0])

%%
% Judits fit to data 
eta = logspace(-2,1,1000);
a = 0.18435;
b = 0.30802;
c = 0.02164;
fR = @(eta) a*eta./(b+eta)+c;
plot(1./eta,fR(eta),'*')

%% Look at gradient of eta
syms Lx Lz x z eta1 eta0 eta(x,z,Lx,Lz,eta0,eta1)
inp = [x,z,Lx,Lz,eta0,eta1];

eta(inp) = eta0 + eta1*exp(-(x/Lx)^2 -(z/Lz)^2);
gradeta = gradient(eta,x);

mf_eta = matlabFunction(eta,inp);
mf_gradeta = matlabFunction(gradeta,inp);


Lx = linspace(0.1,40,10);
Lz = 1;
x = linspace(-40,40,100);
z = 0;
eta0 = 0.0015;
eta1 = 2.4;
[X,LX] = meshgrid(x,Lx);

ETA = eta(X,z,LX,Lz,eta0,eta1);
GRADETA = gradeta(X,z,LX,Lz,eta0,eta1);

nrows = 2;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,X',ETA')

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,X',GRADETA')

%% Fit using symvar
% Data from simulation
Lz = 1;
Lx = [0.5,1,2,5,6,7,8,10,15,20,40];
R = [0.22602171,0.24566273,0.25104229,0.20345172,0.18692082,0.17226913,0.15952006,0.14030094,0.09496913,0.05082819,0.03705537];
% Add zero value
%Lx = [0 Lx];
%R = [0.03 R];

syms x a b c d f(x,a,b,c,d)
inp = [x,a,b,c,d];

% b = 8; % b
% a = 2.4*8; % a
% c = 0.18; % c
% E0 = 0.03;
% f = @(Lx) a/b/b*(c+1)*(Lx/b).^(c).*exp(-(Lx/b).^(c+1)) + E0;
    
f(inp) = a*exp(-(x/b)^(c+1)) + d*x;
gradf(inp) = gradient(f,x);

mf = matlabFunction(f);
mf_grad = matlabFunction(gradf);

params0 = [-2.4,8,0.18,0.03];
cf = @(params) costfunction(params,Lx,R,mf_grad);
params = fminsearch(cf,params0);

Rfit0 = mf_grad(Lx,params0(1),params0(2),params0(3),params0(4));
Rfit = mf_grad(Lx,params(1),params(2),params(3),params(4));

nrows = 1;
ncols = 1;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
plot(hca,Lx,R,'*',Lx,Rfit0,'--',Lx,Rfit,'-')
legend({'data','initial guess','final fit'})

%hca = subplot(nrows,ncols,isub); isub = isub + 1;
%plot(hca,X',GRADETA')


function cost = costfunction(params,dep_var,data,fit_function)
  % cost = costfunction(params,dep_var,fit_data,fit_function)
  % params = [x,z,Lx,Lz,eta0,eta1]
  % dep_Var = Lx
  % fit_data = ER
  % fit_function = matlab function to fit

  a = params(1);
  b = params(2);
  c = params(3);
  d = params(4);
  data_fit = fit_function(dep_var,a,b,c,d);


  colors = mms_colors('1234');
  cost = sum((data(:)-data_fit(:)).^2);
  
  if 0 % plot 
    hca = subplot(1,1,1);
    plot(hca,dep_var,data,dep_var,data_fit)  
    pause(0.1)
  end
end










