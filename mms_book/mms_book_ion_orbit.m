no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

%%
twpe = 20000;
%A0 = no02m.twpelim(twpe).A;
%[X0,Z0] = ndgrid(no02m.xi,no02m.zi);

% Integrate trajectories
pic = no02m.twpelim(twpe).xlim([90 105]).zlim([-5 5]);
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



x0 = [100.6];
y0 = 0;
z0 = [-4];

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


fprintf('[%7.4f, %7.4f, %7.4f]\n',vx0,vy0,vz0)
tstart = 0;
tstop = 200;
m = 1;
q = 1;


        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);

options = odeset('Events',@exitBox,'Events', @(tt,xx) myEventBoxEdge(tt,xx,pic.xi(([1 end]))),'AbsTol',1e-14,'RelTol',1e-7);
EoM = @(t,xyz_vxvyvz) eom_pic_external(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,m,q); 

clear pall
for ipart = 1:numel(x0)  
  r0 = [x0(ipart) 0 z0(ipart)];
  v0 = [vx0(ipart) vy0(ipart) vz0(ipart)];
  x_init = [r0, v0];

  [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);
    
  pall(ipart).t = t;
  pall(ipart).x = x_sol_tmp(:,1);
  pall(ipart).y = x_sol_tmp(:,2);
  pall(ipart).z = x_sol_tmp(:,3);
  pall(ipart).vx = x_sol_tmp(:,4);
  pall(ipart).vy = x_sol_tmp(:,5);
  pall(ipart).vz = x_sol_tmp(:,6);
  pall(ipart).Ex = fEx(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall(ipart).Ey = fEy(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall(ipart).Ez = fEz(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall(ipart).Bx = fBx(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall(ipart).By = fBy(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall(ipart).Bz = fBz(x_sol_tmp(:,1),x_sol_tmp(:,3));
end

% Plot results
ipart = 1;
p = pall(ipart);


Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;

h = setup_subplots(2,2);
isub = 1;

if 0 % Ez 
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,fEz(X,Z))
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,p.x,p.z,'Color','k','linewidth',1)
  hold(hca,'off')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'E_z';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hold(hca,'on')  
  contour(hca,X0,Z0,A0,Alevs,'k')
  hold(hca,'off')
end
if 0 % Ex 
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,fEx(X,Z))
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,p.x,p.z,'Color','k','linewidth',1)
  hold(hca,'off')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'E_x';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hold(hca,'on')  
  contour(hca,X0,Z0,A0,Alevs,'k')
  hold(hca,'off')
  hca.CLim = [-0.4 0.4];
end

if 1 % Ex, with quivers
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,fEx(X,Z))
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,p.x,p.z,'Color','k','linewidth',3)
  hold(hca,'off')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'E_x';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hold(hca,'on')  
  contour(hca,X0,Z0,A0,Alevs,'k')
  hold(hca,'off')

  hold(hca,'on') 
  xx = p.x;
  zz = p.z;  
  quiver(hca,xx(1:quiv_step:end),zz(1:quiv_step:end),quiv_scal*p.Ex(1:quiv_step:end),quiv_scal*p.Ez(1:quiv_step:end),0,'color','k','LineWidth',1)  
  axis(hca,'equal')
  hca.XLim = pic.xi([1 end]);
  hold(hca,'off')
  hca.Box = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'z/d_i';

  irf_legend(hca,sprintf('N_{smooth} = %g', ns*2+1),[0.98 1.04],'k')
end

if 1 % By
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,fBy(X,Z))
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,p.x,p.z,'Color','k','linewidth',3)
  hold(hca,'off')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'B_y';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hold(hca,'on')  
  contour(hca,X0,Z0,A0,Alevs,'k')
  hold(hca,'off')

  hold(hca,'on') 
  xlim = hca.XLim;
  zlim = hca.YLim;
  vxB_x = -p.vz.*p.By + 0*p.vy.*p.Bz;
  vxB_z = p.vx.*p.By - 0*p.vy.*p.Bx;
  xx = p.x;
  zz = p.z;  
  quiver(hca,xx(1:quiv_step:end),zz(1:quiv_step:end),quiv_scal*vxB_x(1:quiv_step:end),quiv_scal*vxB_z(1:quiv_step:end),0,'color','k','LineWidth',1)
  
  vxB_x = 0*-p.vz.*p.By + 1*p.vy.*p.Bz;
  vxB_z = 0*p.vx.*p.By - 1*p.vy.*p.Bx;
  xx = p.x;
  zz = p.z;    
  quiver(hca,xx(1:quiv_step:end),zz(1:quiv_step:end),quiv_scal*vxB_x(1:quiv_step:end),quiv_scal*vxB_z(1:quiv_step:end),0,'color',[0.5 0.5 0.5],'LineWidth',1)
  axis(hca,'equal')
  hca.XLim = pic.xi([1 end]);
  hca.XLim = xlim;
  hca.YLim = zlim;
  hold(hca,'off')


  hca.Box = 'on';
  hca.Layer = 'top';

  hca.Title.String = {sprintf('r_0 = [%7.4f, %7.4f, %7.4f]\n',x0,y0,z0);
                      sprintf('v_0 = [%7.4f, %7.4f, %7.4f]\n',vx0,vy0,vz0)};
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'z/d_i';

  irf_legend(hca,sprintf('N_{smooth} = %g', ns*2+1),[0.98 1.04],'k')
end
if 1 % 3D trajectory
  hca = h(isub); isub = isub + 1;
  % pcolor(hca,X,Z,fBy(X,Z))
  % shading(hca,'flat')
  % hold(hca,'on')
  % plot(hca,p.x,p.z,'Color','k','linewidth',1)
  % hold(hca,'off')
  % hcb = colorbar(hca);
  % hcb.YLabel.String = 'B_y';
  % hca.CLim = max(abs(hca.CLim))*[-1 1];
  % hold(hca,'on')  
  % contour(hca,X0,Z0,A0,Alevs,'k')
  % hold(hca,'off')

  plot3(hca,p.x,p.y,p.z,'Color','k','linewidth',2)

  hold(hca,'on') 
  vxB_x = -p.vz.*p.By + p.vy.*p.Bz;
  vxB_y = p.vz.*p.Bx - p.vx.*p.Bz;
  vxB_z = p.vx.*p.By - p.vy.*p.Bx;
  xx = p.x;
  yy = p.y;
  zz = p.z;  
  % vxB
  quiver3(hca,xx(1:quiv_step:end),yy(1:quiv_step:end),zz(1:quiv_step:end),...
    quiv_scal*vxB_x(1:quiv_step:end),quiv_scal*vxB_y(1:quiv_step:end),quiv_scal*vxB_z(1:quiv_step:end),0,'color','b','LineWidth',1)
  
  % E
  quiver3(hca,xx(1:quiv_step:end),yy(1:quiv_step:end),zz(1:quiv_step:end),...
    1*quiv_scal*p.Ex(1:quiv_step:end),1*quiv_scal*p.Ey(1:quiv_step:end),1*quiv_scal*p.Ez(1:quiv_step:end),0,'color','r','LineWidth',1)
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  view(hca,[0 0 1]);
end
if 0 % L-forces vs time
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.Ex, p.t,-p.vz.*p.By, p.t, p.vy.*p.Bz,'linewidth',1)
  legend(hca,{'E_x','-v_zB_y','v_yB_z'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Cumalative L-forces vs time
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,cumtrapz(p.t,p.Ex), p.t,cumtrapz(p.t,-p.vz.*p.By), p.t, cumtrapz(p.t,p.vy.*p.Bz),'linewidth',1)
  hold(hca,'on')
  plot(hca,p.t,p.vx,'linewidth',1,'color',[0 0 0])
  hold(hca,'off')
  %plot(hca,p.t,cumsum(p.Ex), p.t,cumsum(-p.vz.*p.By), p.t, cumsum(p.vy.*p.Bz),'linewidth',1)
  legend(hca,{'E_x','-v_zB_y','v_yB_z','v_x'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'v_x = (1/m)\int F_x dt';
  hca.XLabel.String = '(t-t_{start})\omega_{ci}';
end
if 0 % position vs time
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.x-p.x(1), p.t,p.z,'linewidth',1)
  legend(hca,{'x-x_{start}','z'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Cumalative L-forces vs x
  hca = h(isub); isub = isub + 1;
  plot(hca,p.z,cumsum(p.Ex), p.z,cumsum(-p.vz.*p.By), p.z, cumsum(p.vy.*p.Bz),'linewidth',1)
  legend(hca,{'E_x','-v_zB_y','v_yB_z'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Ex, with scatter of Ex for bugcheck
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,fEy(X,Z))
  shading(hca,'flat')
  hold(hca,'on')
  scatter(hca,p.x,p.z,5,p.Ey)
  hold(hca,'off')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'E_y';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end


colormap(pic_colors('blue_red'))
c_eval('h(?).FontSize = 14;',1:numel(h))

















