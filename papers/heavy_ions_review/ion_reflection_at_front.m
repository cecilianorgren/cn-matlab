%% ion and electron orbit from the same place...
twpe = 22000;
%A0 = no02m.twpelim(twpe).A;
%[X0,Z0] = ndgrid(no02m.xi,no02m.zi);

% Integrate trajectories
pic = no02m.twpelim(twpe).xlim([20 105]).zlim([-5 5]);
pic = no02m.twpelim(twpe);%.xlim([00 200]);
x = pic.xi;
z = pic.zi;
[X,Z] = ndgrid(x,z);


A0 = pic.twpelim(twpe).A;
[X0,Z0] = ndgrid(pic.xi,pic.zi);

ns = 10;
fEx = griddedInterpolant(X,Z,smooth2(pic.Ex,ns));
fEy = griddedInterpolant(X,Z,smooth2(pic.Ey,ns));
fEz = griddedInterpolant(X,Z,smooth2(pic.Ez,ns));
fBx = griddedInterpolant(X,Z,pic.Bx);
fBy = griddedInterpolant(X,Z,pic.By);
fBz = griddedInterpolant(X,Z,pic.Bz);
fvx = griddedInterpolant(X,Z,pic.vix);
fvy = griddedInterpolant(X,Z,pic.viy);
fvz = griddedInterpolant(X,Z,pic.viz);



x0 = [70];
y0 = 0;
z0 = [5];

ip = 1;

vx0 = fvx(x0',z0');
vy0 = fvy(x0',z0');
vz0 = fvx(x0',z0');


Babs = sqrt(fBx(x0',z0').^2 + fBy(x0',z0').^2 + fBz(x0',z0').^2);
vx0 = (fEy(x0',z0')*fBx(x0',z0')-fEx(x0',z0')*fBy(x0',z0'))/Babs;
vy0 = (fEz(x0',z0')*fBx(x0',z0')-fEx(x0',z0')*fBz(x0',z0'))/Babs;
vz0 = (fEx(x0',z0')*fBy(x0',z0')-fEy(x0',z0')*fBx(x0',z0'))/Babs;



%vx0 = 0;
%vx0 = .0;
%vz0 = -0.2;
%vy0 = 0;

vx0i = vx0 + 1*0.3;
vy0i = vy0;
vz0i = vz0*0;

vx0e = vx0 + 2;
vy0e = vy0 + 0.5;
vz0e = vz0 - 0.5;

vx0e = vx0 + 2;
vy0e = vy0 + 0.5;
vz0e = vz0 - 0.5; % me = 1/400

vx0e = vx0 + 2;
vy0e = vy0 + 0.5;
vz0e = vz0 - 0.5;

fprintf('[%7.4f, %7.4f, %7.4f]\n',vx0,vy0,vz0)
tstart = 0;
tstop = 100;
mi = 1;
qi = 1;
me = 1/1836;
qe = -1;


        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);

%options = odeset('Events',@exitBox,'Events', @(tt,xx) myEventBoxEdge(tt,xx,pic.xi(([1 end]))),'AbsTol',1e-14,'RelTol',1e-7);
options = odeset('AbsTol',1e-7,'RelTol',1e-7);
EoMi = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,mi,qi); 
EoMe = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,me,qe); 

clear pall_e pall_i
tic
for ipart = 1:numel(x0)  
  r0 = [x0(ipart) 0 z0(ipart)];
  v0i = [vx0i(ipart) vy0i(ipart) vz0i(ipart)];
  v0e = [vx0e(ipart) vy0e(ipart) vz0e(ipart)];
  x_init_i = [r0, v0i];
  x_init_e = [r0, v0e];

  [ti,x_sol_tmpi] = ode45(EoMi,[tstart tstop],x_init_i,options);  
    
  pall_i(ipart).t = ti;
  pall_i(ipart).x = x_sol_tmpi(:,1);
  pall_i(ipart).y = x_sol_tmpi(:,2);
  pall_i(ipart).z = x_sol_tmpi(:,3);
  pall_i(ipart).vx = x_sol_tmpi(:,4);
  pall_i(ipart).vy = x_sol_tmpi(:,5);
  pall_i(ipart).vz = x_sol_tmpi(:,6);
  pall_i(ipart).Ex = fEx(x_sol_tmpi(:,1),x_sol_tmpi(:,3));
  pall_i(ipart).Ey = fEy(x_sol_tmpi(:,1),x_sol_tmpi(:,3));
  pall_i(ipart).Ez = fEz(x_sol_tmpi(:,1),x_sol_tmpi(:,3));
  pall_i(ipart).Bx = fBx(x_sol_tmpi(:,1),x_sol_tmpi(:,3));
  pall_i(ipart).By = fBy(x_sol_tmpi(:,1),x_sol_tmpi(:,3));
  pall_i(ipart).Bz = fBz(x_sol_tmpi(:,1),x_sol_tmpi(:,3));

  if 0
  [te,x_sol_tmpe] = ode45(EoMe,[tstart tstop],x_init_e,options);
  pall_e(ipart).t = te;
  pall_e(ipart).x = x_sol_tmpe(:,1);
  pall_e(ipart).y = x_sol_tmpe(:,2);
  pall_e(ipart).z = x_sol_tmpe(:,3);
  pall_e(ipart).vx = x_sol_tmpe(:,4);
  pall_e(ipart).vy = x_sol_tmpe(:,5);
  pall_e(ipart).vz = x_sol_tmpe(:,6);
  pall_e(ipart).Ex = fEx(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Ey = fEy(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Ez = fEz(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Bx = fBx(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).By = fBy(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Bz = fBz(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  end
end
toc

% Plot results
ipart = 1;
%p_e = pall_e(ipart);
p_i = pall_i(ipart);


Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;

h = setup_subplots(2,1);
isub = 1;
tlim = [0 Inf];

if 1 % 2D trajectory, xz
  hca = h(isub); isub = isub + 1;
  

  Ey = pic.Ey;
  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  pcolor(hca,xi,zi,smooth2(Ey,3)')
  shading(hca,'flat')
  
  hold(hca,'on')
  Ay_levels = min(Ay(:)):1:max(Ay(:));

  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2))); 
  plot(hca,p_i.x(itpi),p_i.z(itpi),'linewidth',2)  
  hold(hca,'off') 
  
  %axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [40 160];
end
if 1 % 2D trajectory, xy
  hca = h(isub); isub = isub + 1;

  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2))); 
  plot(hca,p_i.x(itpi),p_i.y(itpi),'linewidth',2)  

  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [40 160];
end
if 0 % 3D trajectory
  hca = h(isub); isub = isub + 1;

  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  Ay_levels = min(Ay(:)):1:max(Ay(:));
  %contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
  %contour3(hca,Xi,Yi,Zi,Ay,Ay_levels,'color','k') 
  %%
  

  plot3(hca,p_i.x,p_i.y,p_i.z,p_e.x,p_e.y,p_e.z,'linewidth',2)  
  %plot3(hca,'Color','b','linewidth',2)
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
end



colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 14;',1:numel(h))
































%% ion and electron orbit from the same place... including time
twpe = 18000;
%A0 = no02m.twpelim(twpe).A;
%[X0,Z0] = ndgrid(no02m.xi,no02m.zi);

% Integrate trajectories
%pic = no02m.twpelim(twpe).xlim([20 105]).zlim([-5 5]);
pic = no02m.twpelim([14000 25000]);%.xlim([100 140]);
x = pic.xi;
z = pic.zi;
t = pic.twci;
[X,Z,T] = ndgrid(x,z,t);


%A0 = pic.twpelim(twpe).A;
%[X0,Z0] = ndgrid(pic.xi,pic.zi);
tic
ns = 10;
Ex = pic.Ex; for it = 1:pic.nt; Ex(:,:,it) = smooth2(Ex(:,:,it),ns); end
Ey = pic.Ey; for it = 1:pic.nt; Ey(:,:,it) = smooth2(Ey(:,:,it),ns); end
Ez = pic.Ez; for it = 1:pic.nt; Ez(:,:,it) = smooth2(Ez(:,:,it),ns); end
fEx = griddedInterpolant(X,Z,T,Ex);
fEy = griddedInterpolant(X,Z,T,Ey);
fEz = griddedInterpolant(X,Z,T,Ez);
fBx = griddedInterpolant(X,Z,T,pic.Bx);
fBy = griddedInterpolant(X,Z,T,pic.By);
fBz = griddedInterpolant(X,Z,T,pic.Bz);
fvx = griddedInterpolant(X,Z,T,pic.vix);
fvy = griddedInterpolant(X,Z,T,pic.viy);
fvz = griddedInterpolant(X,Z,T,pic.viz);
toc

%%

x0 = [85];
y0 = 0;
z0 = [2];

tstart = 60;
tstop = 125;
mi = 1;
qi = 1;
me = 1/100;
qe = -1;

ip = 1;

vx0 = fvx(x0',z0',tstart);
vy0 = fvy(x0',z0',tstart);
vz0 = fvx(x0',z0',tstart);


Babs = sqrt(fBx(x0',z0',tstart).^2 + fBy(x0',z0',tstart).^2 + fBz(x0',z0',tstart).^2);
vx0 = (fEy(x0',z0',tstart)*fBx(x0',z0',tstart)-fEx(x0',z0',tstart)*fBy(x0',z0',tstart))/Babs;
vy0 = (fEz(x0',z0',tstart)*fBx(x0',z0',tstart)-fEx(x0',z0',tstart)*fBz(x0',z0',tstart))/Babs;
vz0 = (fEx(x0',z0',tstart)*fBy(x0',z0',tstart)-fEy(x0',z0',tstart)*fBx(x0',z0',tstart))/Babs;

%vx0 = 0;
%vx0 = .0;
%vz0 = -0.2;
%vy0 = 0;


fprintf('[%7.4f, %7.4f, %7.4f]\n',vx0,vy0,vz0)



        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);

%options = odeset('Events',@exitBox,'Events', @(tt,xx) myEventBoxEdge(tt,xx,pic.xi(([1 end]))),'AbsTol',1e-14,'RelTol',1e-7);
options = odeset('AbsTol',1e-14,'RelTol',1e-7);
EoMi = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant_t(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,mi,qi); 
EoMe = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant_t(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,me,qe); 

clear pall_e pall_i
tic
for ipart = 1:numel(x0)  
  r0 = [x0(ipart) 0 z0(ipart)];
  v0 = [vx0(ipart) vy0(ipart) vz0(ipart)];
  x_init = [r0, v0];

  [ti,x_sol_tmpi] = ode45(EoMi,[tstart tstop],x_init,options);
  
    
  pall_i(ipart).t = ti;
  pall_i(ipart).x = x_sol_tmpi(:,1);
  pall_i(ipart).y = x_sol_tmpi(:,2);
  pall_i(ipart).z = x_sol_tmpi(:,3);
  pall_i(ipart).vx = x_sol_tmpi(:,4);
  pall_i(ipart).vy = x_sol_tmpi(:,5);
  pall_i(ipart).vz = x_sol_tmpi(:,6);
  pall_i(ipart).Ex = fEx(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);
  pall_i(ipart).Ey = fEy(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);
  pall_i(ipart).Ez = fEz(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);
  pall_i(ipart).Bx = fBx(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);
  pall_i(ipart).By = fBy(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);
  pall_i(ipart).Bz = fBz(x_sol_tmpi(:,1),x_sol_tmpi(:,3),ti);

  if 0
    [te,x_sol_tmpe] = ode45(EoMe,[tstart tstop],x_init,options);
    pall_e(ipart).t = te;
    pall_e(ipart).x = x_sol_tmpe(:,1);
    pall_e(ipart).y = x_sol_tmpe(:,2);
    pall_e(ipart).z = x_sol_tmpe(:,3);
    pall_e(ipart).vx = x_sol_tmpe(:,4);
    pall_e(ipart).vy = x_sol_tmpe(:,5);
    pall_e(ipart).vz = x_sol_tmpe(:,6);
    pall_e(ipart).Ex = fEx(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
    pall_e(ipart).Ey = fEy(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
    pall_e(ipart).Ez = fEz(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
    pall_e(ipart).Bx = fBx(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
    pall_e(ipart).By = fBy(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
    pall_e(ipart).Bz = fBz(x_sol_tmpe(:,1),x_sol_tmpe(:,3),te);
  end
end
toc

%% Plot results
ipart = 1;
%p_e = pall_e(ipart);
p_i = pall_i(ipart);

fontsize = 16; 
Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;

h = setup_subplots(4,2);
isub = 1;
xlim = [40 160];

if 1 % 2D trajectory for a few different times
  dt = 10;
  times = tstart:dt:tstop;  
  Ay_all = pic.twcilim(times,'exact').A;
  Ay_levels = min(Ay_all(:)):1:max(Ay_all(:));
  for it = 1:numel(times)
    hca = h(isub); isub = isub + 1;

    time = times(it);
    Ay = Ay_all(:,:,it);
    xi = pic.xi;
    zi = pic.zi;
    [Xi,Zi] = ndgrid(xi,zi);
    contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
    
    hold(hca,'on')
    itpi = intersect(find(p_i.t>time-0.5*dt),find(p_i.t<time+0.5*dt));
    plot(hca,p_i.x(itpi),p_i.z(itpi),'linewidth',2)  
    hold(hca,'off') 
    
    hold(hca,'off')
    axis(hca,'equal')
    hca.Box = 'on';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'y/d_i';
    hca.ZLabel.String = 'z/d_i';
    %view(hca,[0 1 0]);
    hca.XLim = xlim;
    irf_legend(hca,{sprintf('tw_{ci} = %g',time)},[0.02 0.98],'color','k','fontsize',fontsize)
  end
end
if 0 % 2D trajectory
  hca = h(isub); isub = isub + 1;

  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  Ay_levels = min(Ay(:)):1:max(Ay(:));
  contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
  %contour3(hca,Xi,Yi,Zi,Ay,Ay_levels,'color','k') 
  %%
  hold(hca,'on')

  %plot3(hca,p_i.x,p_i.y,p_i.z,p_e.x,p_e.y,p_e.z,'linewidth',2)  
  %plot3(hca,'Color','b','linewidth',2)
  plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
end
if 0 % 2D trajectory
  hca = h(isub); isub = isub + 1;

  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  Ay_levels = min(Ay(:)):1:max(Ay(:));
  contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
  %contour3(hca,Xi,Yi,Zi,Ay,Ay_levels,'color','k') 
  %%
  hold(hca,'on')

  %plot3(hca,p_i.x,p_i.y,p_i.z,p_e.x,p_e.y,p_e.z,'linewidth',2)  
  %plot3(hca,'Color','b','linewidth',2)
  plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
end
if 0 % 3D trajectory
  hca = h(isub); isub = isub + 1;

  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  Ay_levels = min(Ay(:)):1:max(Ay(:));
  %contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
  %contour3(hca,Xi,Yi,Zi,Ay,Ay_levels,'color','k') 
  %%
  

  plot3(hca,p_i.x,p_i.y,p_i.z,p_e.x,p_e.y,p_e.z,'linewidth',2)  
  %plot3(hca,'Color','b','linewidth',2)
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
end


colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 14;',1:numel(h))


















