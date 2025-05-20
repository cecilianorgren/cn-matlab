no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%no02m = PIC('/Users/cecilianorgren/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

%%
twpe = 20000;
%A0 = no02m.twpelim(twpe).A;
%[X0,Z0] = ndgrid(no02m.xi,no02m.zi);

% Integrate trajectories
pic = no02m.twpelim(twpe).xlim([80 105]).zlim([-5 5]);
x = pic.xi;
z = pic.zi;
[X,Z] = ndgrid(x,z);


intEz = pic.intEz;
phi = -intEz;
A0 = pic.A;
AX = pic.Ax;
[X0,Z0] = ndgrid(pic.xi,pic.zi);



ns = 10;
fPhi = griddedInterpolant(X,Z,smooth2(phi,ns));
fAx = griddedInterpolant(X,Z,smooth2(AX,ns));
fAy = griddedInterpolant(X,Z,smooth2(A0,ns));
fEx = griddedInterpolant(X,Z,smooth2(pic.Ex,ns));
fEy = griddedInterpolant(X,Z,smooth2(pic.Ey,ns));
fEz = griddedInterpolant(X,Z,smooth2(pic.Ez,ns));
fBx = griddedInterpolant(X,Z,pic.Bx);
fBy = griddedInterpolant(X,Z,pic.By);
fBz = griddedInterpolant(X,Z,pic.Bz);
fvx = griddedInterpolant(X,Z,pic.vix);
fvy = griddedInterpolant(X,Z,pic.viy);
fvz = griddedInterpolant(X,Z,pic.viz);
%%


x0 = [99.0];
y0 = 0;
z0 = [-3];

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
tstop = 20;
mi = 1;
qe = 1;


        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);

%options = odeset('Events',@exitBox,'Events', @(tt,xx) myEventBoxEdge(tt,xx,pic.xi(([1 end]))),'AbsTol',1e-14,'RelTol',1e-7);
options = odeset('AbsTol',1e-14,'RelTol',1e-7);
EoM = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,mi,qe); 

clear pall_e
for ipart = 1:numel(x0)  
  r0 = [x0(ipart) 0 z0(ipart)];
  v0 = [vx0(ipart) vy0(ipart) vz0(ipart)];
  x_init = [r0, v0];

  [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);
    
  pall_e(ipart).t = t;
  pall_e(ipart).x = x_sol_tmp(:,1);
  pall_e(ipart).y = x_sol_tmp(:,2);
  pall_e(ipart).z = x_sol_tmp(:,3);
  pall_e(ipart).vx = x_sol_tmp(:,4);
  pall_e(ipart).vy = x_sol_tmp(:,5);
  pall_e(ipart).vz = x_sol_tmp(:,6);
  pall_e(ipart).Phi = fPhi(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Ax = fAx(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Ay = fAy(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Ex = fEx(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Ey = fEy(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Ez = fEz(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Bx = fBx(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).By = fBy(x_sol_tmp(:,1),x_sol_tmp(:,3));
  pall_e(ipart).Bz = fBz(x_sol_tmp(:,1),x_sol_tmp(:,3));
end

% Plot results
ipart = 1;
p = pall_e(ipart);


Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;
vlim = 2;
fontsize = 16;

h = setup_subplots(3,2);
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

  irf_legend(hca,sprintf('N_{smooth} = %g', ns*2+1),[0.98 1.06],'k')
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

  irf_legend(hca,sprintf('N_{smooth} = %g', ns*2+1),[0.98 1.06],'k')
end
if 0 % 3D trajectory
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
if 1 % L-forces vs time
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.Ex, p.t,-p.vz.*p.By, p.t, p.vy.*p.Bz,'linewidth',1)
  %legend(hca,{'E_x','-v_zB_y','v_yB_z'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'F_x';
  hca.XLabel.String = 't\omega_{ci}';
  irf_legend(hca,{'E_x','-v_zB_y','v_yB_z'},[0.02 0.15],'fontsize',fontsize)
end
if 0 % Cumalative L-forces vs time
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
  hca.XLabel.String = 't ()';
end
if 1 % Cumalative L-forces vs x
  hca = h(isub); isub = isub + 1;
  plot(hca,p.x,cumtrapz(p.t,p.Ex), ...
    p.x,cumtrapz(p.t,-p.vz.*p.By), ...
    p.x, cumtrapz(p.t,p.vy.*p.Bz),'linewidth',1)
  hold(hca,'on')
  %plot(hca,p.x,p.vx,'linewidth',1,'color',[0 0 0])
  hold(hca,'off')
  %plot(hca,p.t,cumsum(p.Ex), p.t,cumsum(-p.vz.*p.By), p.t, cumsum(p.vy.*p.Bz),'linewidth',1)
  %legend(hca,{'E_x','-v_zB_y','v_yB_z','v_x'},'Box','off','Location','best')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'v_x = (1/m)\int F_x dt';
  hca.XLabel.String = 'x (d_i)';
  irf_legend(hca,{'E_x','-v_zB_y','v_yB_z'},[0.98 0.15],'fontsize',fontsize)
end
if 0 % Cumalative L-forces vs x, one more integration to get x
  hca = h(isub); isub = isub + 1;
  x0 = p.x(1);
  int1 = cumtrapz(p.t,cumtrapz(p.t,p.Ex));
  int2 = cumtrapz(p.t,cumtrapz(p.t,-p.vz.*p.By));
  int3 = cumtrapz(p.t,cumtrapz(p.t,p.vy.*p.Bz));
  plot(hca,p.x-x0, int1, p.x-x0, int2, p.x-x0, int3,'linewidth',1)
  hold(hca,'on')
  plot(hca,p.x-x0, int1 + int2 + int3,'linewidth',1,'color',[0 0 0])
  hold(hca,'off')
  %plot(hca,p.t,cumsum(p.Ex), p.t,cumsum(-p.vz.*p.By), p.t, cumsum(p.vy.*p.Bz),'linewidth',1)
  legend(hca,{'E_x','-v_zB_y','v_yB_z','v_x'},'Box','off','Location','best')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLabel.String = 'l_x = (1/m)\int \int F_x dt';
  hca.XLabel.String = 'x-x_(t=0) (d_i)';
  axis(hca,'equal')
  axis(hca,'square')
end
if 1 % vx, vz
  hca = h(isub); isub = isub + 1;
  plot(hca,p.vx,p.vz,p.vx(1),p.vz(1),'go','linewidth',2)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  axis(hca,'equal')
  hca.XLim = vlim*[-1 1];
  hca.YLim = vlim*[-1 1];
  hca.XTick = hca.YTick;
end
if 0 % vx, vy
  hca = h(isub); isub = isub + 1;
  plot(hca,p.vx,p.vy,p.vx(1),p.vy(1),'go','linewidth',2)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  axis(hca,'equal')
  hca.XLim = vlim*[-1 1];
  hca.YLim = vlim*[-1 1];
  hca.XTick = hca.YTick;
end
if 0 % Ax, vx
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.vx,p.t,p.Ax,p.t,p.vx+p.Ax,'linewidth',1)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'vx, Ax';  
  legend(hca,{'v_x','A_x','v_x+A_x'},'Box','off','location','best')  
end
if 0 % Ay, vy
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.vy,p.t,p.Ay,p.t,p.vy+p.Ay,'linewidth',1)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'vy, Ay';  
  legend(hca,{'v_y','A_y','v_y+A_y'},'Box','off','location','best')  
end
if 0 % position vs time
  hca = h(isub); isub = isub + 1;
  plot(hca,p.t,p.x-p.x(1), p.t,p.z,'linewidth',1)
  legend(hca,{'x-x_{start}','z'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % injection model, v^2
  hca = h(isub); isub = isub + 1;
  v2 = sqrt(p.vx.^2 + p.vy.^2 + p.vz.^2);
  plot(hca,p.t,v2,'linewidth',1)
  %legend(hca,{'v_x^{mod}','v_x'},'Box','off')  
  hca.XLabel.String = 't';
  hca.YLabel.String = 'v^2';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

if 0 % injection model, phi, Ax
  hca = h(isub); isub = isub + 1;
  v2 = sqrt(p.vx.^2 + 0*p.vy.^2 + p.vz.^2);
  plot(hca,p.t,p.Ax, p.t,p.Phi, p.t,-v2/2,'linewidth',1)
  legend(hca,{'A_x','\phi','-v^2/2'},'Box','off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % injection model, Phi/Ax
  hca = h(isub); isub = isub + 1;
  toplot = p.Phi./p.Ax;
  toplot(abs(toplot)>prctile(abs(toplot),96)) = NaN;
  plot(hca,p.t,toplot,'linewidth',1)  
  hca.XLabel.String = 't';
  hca.YLabel.String = '\phi/A_x';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % injection model
  hca = h(isub); isub = isub + 1;
  v2 = sqrt(p.vx.^2 + p.vy.^2 + p.vz.^2);
  toplot = -(p.Ax./p.Phi).*v2/2;
  toplot(abs(toplot)>prctile(abs(toplot),96)) = NaN;
  plot(hca,p.t,toplot, p.t,p.vx,'linewidth',1)
  legend(hca,{'v_x^{mod}','v_x'},'Box','off')  
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
c_eval('h(?).FontSize = 16;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); 
c_eval('hl(?).LineWidth = 2;',1:numel(hl))



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



x0 = [98];
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

vx0i = vx0 + 0*0.2;
vy0i = vy0;
vz0i = vz0;

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
tstop = 60;
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
  [te,x_sol_tmpe] = ode45(EoMe,[tstart tstop],x_init_e,options);
    
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
toc

%% Plot results
ipart = 1;
p_e = pall_e(ipart);
p_i = pall_i(ipart);


Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;

h = setup_subplots(2,1);
isub = 1;
tlim = [0 40];

if 1 % 2D trajectory
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
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  plot(hca,p_i.x(itpi),p_i.z(itpi),p_e.x(itpe),p_e.z(itpe),'linewidth',2)  
  
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [40 160];
end
if 1 % 3D trajectory
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

x0 = [102];
y0 = 0;
z0 = [-4];

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
  [te,x_sol_tmpe] = ode45(EoMe,[tstart tstop],x_init,options);
    
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
toc

% Plot results
ipart = 1;
p_e = pall_e(ipart);
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
  Ay_all = pic.A;
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
    itpe = intersect(find(p_e.t>time-0.5*dt),find(p_e.t<time+0.5*dt));    
    plot(hca,p_i.x(itpi),p_i.z(itpi),p_e.x(itpe),p_e.z(itpe),'linewidth',2)  
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






























%% 2 ions, with and without the Hall fields
twpe = 20000;
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



x0 = [99];
y0 = 0;
z0 = [-6];

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

vx0i = vx0;
%vx0i = 0;
vy0i = vy0;
vz0i = vz0;

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

fEz_nohall = @(x,z)0*x;
fEx_nohall = @(x,z)0*x;
fBy_nohall = @(x,z)0*x;
        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);

%options = odeset('Events',@exitBox,'Events', @(tt,xx) myEventBoxEdge(tt,xx,pic.xi(([1 end]))),'AbsTol',1e-14,'RelTol',1e-7);
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
EoMi = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant(t,xyz_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,mi,qi); 
EoMi_noHall = @(t,xyz_vxvyvz) eom_pic_external_griddedInterpolant(t,xyz_vxvyvz,fEx_nohall,fEy,fEz_nohall,fBx,fBy_nohall,fBz,mi,qi); 

clear pall_e pall_i
tic
for ipart = 1:numel(x0)  
  r0 = [x0(ipart) 0 z0(ipart)];
  v0i = [vx0i(ipart) vy0i(ipart) vz0i(ipart)];
  v0e = [vx0i(ipart) vy0i(ipart) vz0i(ipart)];
  x_init_i = [r0, v0i];
  x_init_e = [r0, v0e];

  [ti,x_sol_tmpi] = ode45(EoMi,[tstart tstop],x_init_i,options);
  [te,x_sol_tmpe] = ode45(EoMi_noHall,[tstart tstop],x_init_e,options);
    
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

  pall_e(ipart).t = te;
  pall_e(ipart).x = x_sol_tmpe(:,1);
  pall_e(ipart).y = x_sol_tmpe(:,2);
  pall_e(ipart).z = x_sol_tmpe(:,3);
  pall_e(ipart).vx = x_sol_tmpe(:,4);
  pall_e(ipart).vy = x_sol_tmpe(:,5);
  pall_e(ipart).vz = x_sol_tmpe(:,6);
  pall_e(ipart).Ex = fEx_nohall(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Ey = fEy(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Ez = fEz_nohall(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Bx = fBx(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).By = fBy_nohall(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
  pall_e(ipart).Bz = fBz(x_sol_tmpe(:,1),x_sol_tmpe(:,3));
end
toc

%% Plot results
ipart = 1;
p_e = pall_e(ipart);
p_i = pall_i(ipart);


Alevs = min(A0(:)):1:max(A0(:));
quiv_step = 100;
quiv_scal = 6;

fontsize = 16;

h = setup_subplots(3,1);
isub = 1;
tlim = [0 55];

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
  [Xi,Zi] = ndgrid(xi,zi);
  
  Ex = fEx(Xi,Zi);
  pcolor(hca,xi,zi,Ex')
  shading(hca,'flat')
  hca.CLim = 0.25*[-1 1];
  hca.Children = circshift(hca.Children,-1);

  %plot3(hca,p_i.x,p_i.y,p_i.z,p_e.x,p_e.y,p_e.z,'linewidth',2)  
  %plot3(hca,'Color','b','linewidth',2)
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  plot(hca,p_i.x(itpi),p_i.z(itpi),p_e.x(itpe),p_e.z(itpe),'linewidth',2)  
  
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [60 140];
  colormap(hca,pic_colors('blue_red'))
end
if 0 % 2D trajectory, with time-colored plots
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
  
  step = 20;
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  tic
  scatter(hca,p_i.x(itpi(1:step:end)),p_i.z(itpi(1:step:end)),5,p_i.t(itpi(1:step:end)))  
  scatter(hca,p_e.x(itpe(1:step:end)),p_e.z(itpe(1:step:end)),5,p_e.t(itpe(1:step:end)))
  toc
  hcb = colorbar(hca);
  hcb.YLabel.String = 't\omega_{pi}';
  colormap(hca,pic_colors('candy4'))
  
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [60 140];
end
if 0 % 2D trajectory, with time-colored plots of vx
  hca = h(isub); isub = isub + 1;

  Ay = pic.A;
  xi = pic.xi;
  zi = pic.zi;
  [Xi,Yi,Zi] = ndgrid(xi,0,zi);
  Ay_levels = min(Ay(:)):1:max(Ay(:));
  contour(hca,xi,zi,Ay',Ay_levels,'color','k') 
  %contour3(hca,Xi,Yi,Zi,Ay,Ay_levels,'color','k') 
  
  hold(hca,'on')
  
  step = 20;
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  tic
  scatter(hca,p_i.x(itpi(1:step:end)),p_i.z(itpi(1:step:end)),5,abs(p_i.vx(itpi(1:step:end))))  
  scatter(hca,p_e.x(itpe(1:step:end)),p_e.z(itpe(1:step:end)),5,abs(p_e.vx(itpe(1:step:end))))
  toc
  hcb = colorbar(hca);
  hcb.YLabel.String = '|v_x|';
  colormap(hca,pic_colors('candy4'))
  hca.CLim = [0 2];
  
  %plot(hca,p_i.x,p_i.z,p_e.x,p_e.z,'linewidth',2)  

  hold(hca,'off') 
  
  hold(hca,'off')
  axis(hca,'equal')
  hca.Box = 'on';
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = [60 140];
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
if 0 % t(x)
  hca = h(isub); isub = isub + 1;

 
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  plot(hca,p_i.x(itpi),p_i.t(itpi),p_e.x(itpe),p_e.t(itpe),'linewidth',2)  
  
  hca.Box = 'on';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  %view(hca,[0 1 0]);
  hca.YLim = tlim;
  irf_legend(hca,{'Normal Hall fields','No Hall fields'}',[0.02 0.02],'fontsize',fontsize)
end
if 0 % t(vx)
  hca = h(isub); isub = isub + 1;

 
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  plot(hca,p_i.vx(itpi),p_i.t(itpi),p_e.vx(itpe),p_e.t(itpe),'linewidth',2)  
  
  hca.Box = 'on';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'v_x/d_i';
  %view(hca,[0 1 0]);
  hca.YLim = tlim;
  irf_legend(hca,{'Normal Hall fields','No Hall fields'}',[0.02 0.02],'fontsize',fontsize)
end
if 1 % x(t)
  hca = h(isub); isub = isub + 1;

 
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));    
  plot(hca,p_i.t(itpi),p_i.x(itpi),p_e.t(itpe),p_e.x(itpe),'linewidth',2)  
  
  hca.Box = 'on';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = tlim;
end
if 1 % Fx(t)
  hca = h(isub); isub = isub + 1;

 
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));  
  Ex1 = p_i.Ex(itpi);
  vzBy1 = -p_i.vz(itpi).*p_i.By(itpi);
  vyBz1 = p_i.vy(itpi).*p_i.Bz(itpi);
  Ex2 = p_e.Ex(itpe);
  vzBy2 = -p_e.vz(itpe).*p_e.By(itpe);
  vyBz2 = p_e.vy(itpe).*p_e.Bz(itpe);
  colors = hca.ColorOrder;
  hca.ColorOrder = colors;
  plot(hca,p_i.t(itpi),Ex1, p_i.t(itpi),vzBy1, p_i.t(itpi),vyBz1, 'linewidth',2)  
  hold(hca,'on')
  hca.ColorOrder = colors;
  plot(hca,p_e.t(itpe),Ex2, p_e.t(itpe),vzBy2, p_e.t(itpe),vyBz2, 'linewidth',2,'linestyle','--')  
  hold(hca,'off')
  hca.Box = 'on';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = tlim;
  irf_legend(hca,{'E_x','-v_yB_x','v_xB_y'},[0.02 0.98],'fontsize',16)
end
if 1 % W(t)
  hca = h(isub); isub = isub + 1;

 
  
  itpi = intersect(find(p_i.t>tlim(1)),find(p_i.t<tlim(2)));
  itpe = intersect(find(p_e.t>tlim(1)),find(p_e.t<tlim(2)));  
  Ex1 = p_i.Ex(itpi); x1 = p_i.x(itpi); 
  Ey1 = p_i.Ey(itpi); y1 = p_i.y(itpi);
  Ez1 = p_i.Ez(itpi); z1 = p_i.z(itpi);
  Ex2 = p_e.Ex(itpe); x2 = p_e.x(itpe);
  Ey2 = p_e.Ey(itpe); y2 = p_e.y(itpe);
  Ez2 = p_e.Ez(itpe); z2 = p_e.z(itpe);
  Wx1 = cumtrapz(x1,Ex1);
  Wy1 = cumtrapz(y1,Ey1);
  Wz1 = cumtrapz(z1,Ez1);
  Wx2 = cumtrapz(x2,Ex2);
  Wy2 = cumtrapz(y2,Ey2);
  Wz2 = cumtrapz(z2,Ez2);

  colors = hca.ColorOrder;
  hca.ColorOrder = colors;
  plot(hca,p_i.t(itpi),Wx1, p_i.t(itpi),Wy1, p_i.t(itpi),Wz1, 'linewidth',2)  
  hold(hca,'on')
  hca.ColorOrder = colors;
  plot(hca,p_e.t(itpe),Wx2, p_e.t(itpe),Wy2, p_e.t(itpe),Wz2, 'linewidth',2,'linestyle','--')  
  hold(hca,'off')
  hca.Box = 'on';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  %view(hca,[0 1 0]);
  hca.XLim = tlim;
  irf_legend(hca,{'W_x','W_y','W_z'},[0.02 0.98],'fontsize',16)
end

h(1).Title.String = sprintf('x_0 = %g, z_0 = %g',x0,z0);


c_eval('h(?).FontSize = 14;',1:numel(h))
































