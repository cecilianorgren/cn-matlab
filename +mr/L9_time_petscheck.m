  units = irf_units;
Er = 1e-3;
vA = 1000e3;


Bz = @(x,t) (units.c/vA)*Er*mr.L9_time_petscheck_Er(t).*(t+x/vA);
Bz = @(x,t) (units.c/vA)*Er*mr.L9_time_petscheck_Er(t);

t = linspace(0,2,50);
x = linspace(-2*vA,2*vA,100);
[X,T] = meshgrid(x,t);

pcolor(X,T,Bz(X,T-X/vA)); shading flat
%%
plot(x,x*0)
set(gca,'ylim',[min(min(Bz(X,T+X/vA))) max(max(Bz(X,T+X/vA)))])
set(gca,'xlim',[0 2*vA])
for it = 2:numel(t)
  plot(x,Bz(x,t(it)-x/vA))  
  set(gca,'xlim',[0 2*vA])
  pause(0.1)
end

%%
x_pos = linspace(0,2*vA,100);
x_neg = linspace(-2*vA,0,100);
x = [x_neg x_pos];
plot(x,x*0)

Bz_pos = @(x,t) (units.c/vA)*Er*mr.L9_time_petscheck_Er(t);
Bz_neg = @(x,t) -(units.c/vA)*Er*mr.L9_time_petscheck_Er(t);


set(gca,'ylim',[0 1])
set(gca,'xlim',[-2*vA 2*vA]);
for it = 2:numel(t)
  plot([x_neg 0],[Bz_neg(x_neg,t(it)+x_neg/vA) 0],[0 x_pos],[0 Bz_pos(x_pos,t(it)-x_pos/vA)])
  set(gca,'ylim',0.4*[-1 1])
  pause(0.1)
end

%% plot a few timesteps
units = irf_units;
Er = 1e-3;
vA = 1000e3;


t = linspace(0.2,2,4);
t = [0.2 0.6 1 1.4];

x_pos = linspace(0,2*vA,100);
x_neg = linspace(-2*vA,0,100);
x = [x_neg x_pos];

Bz_pos = @(x,t) (units.c/vA)*Er*mr.L9_time_petscheck_Er(t);
Bz_neg = @(x,t) -(units.c/vA)*Er*mr.L9_time_petscheck_Er(t);

colors = mms_colors('matlab');

hca = subplot(1,1,1);
plot(x/vA,x*0); 
hold('on')
for it = 1:numel(t)
  lines = plot([x_neg 0]/vA,[Bz_neg(x_neg,t(it)+x_neg/vA) 0],[0 x_pos]/vA,[0 Bz_pos(x_pos,t(it)-x_pos/vA)]);
  lines(1).Color = colors(it,:);
  lines(2).Color = colors(it,:);
end
hold('off')
set(gca,'ylim',0.4*[-1 1])
set(gca,'xlim',[-2*vA 2*vA]/vA);
hca.XLabel.String = 'x/Tv_A';
hca.YLabel.String = 'B_z';
hca.FontSize = 16;
hca.Title.String = sprintf('E_r = %.0f mV/m, v_A = %.0f km/s',Er*1e3,vA*1e-3);
legend({sprintf('t = %.2f',t(1)),sprintf('t = %.2f',t(2)),sprintf('t = %.2f',t(3)),sprintf('t = %.2f',t(4))},'location','northwest')

%%
t = linspace(0,1.2,100);
hca = subplot(1,1,1);
plot(t,mr.L9_time_petscheck_Er(t),t,cumtrapz(t,mr.L9_time_petscheck_Er(t)));
hca.XLabel.String = 't';
hca.YLabel.String = 'E_r,  \Phi/L_X';
legend('E_r','\Phi/L_X');
hca.FontSize = 16;
%hca.Title.String = 'Reconnected flux';
