%% Speiser motions
clear
close

% Constants
q=1;m=1;
nt=29000;
dt=.002;

% initial conditions
x0=0;y0=5;z0=0;
r0=[x0,y0,z0]; % position vector

vz0s=[0.1 1 2 3 4 9 0.01 .5 1 2 4 10];
vy0s=[0 0 0 0 0 0 0 0 0 0];
vx0s=[0 0 0 0 0 0 0 0 0 0];

% Bx = y, Bz=By=1 
% Electric field, Ex=Ey=0
Ez=-0.1;

% intialise variables
x = zeros(nt,1); x(1) = x0;
y = zeros(nt,1); y(1) = y0;
z = zeros(nt,1); z(1) = z0;

for ii=[6]     
    vx = zeros(nt,1); vx(1) = vx0s(ii);
    vy = zeros(nt,1); vy(1) = vy0s(ii);
    vz = zeros(nt,1); vz(1) = vz0s(ii);
    
    for k=2:(nt)
      % Acceleration, depend on velocity, magnetic and electric field
      ax = 0;
      ay = q/m*vz(k-1)*y(k-1);
      az = q/m*(Ez-vy(k-1)*y(k-1));

      % Velocities
      vx(k) = vx(k-1) + ax*dt;
      vy(k) = vy(k-1) + ay*dt;
      vz(k) = vz(k-1) + az*dt;

      % Positions
      x(k) = x(k-1) + vx(k-1)*dt + .5*ax*dt^2;
      y(k) = y(k-1) + vy(k-1)*dt + .5*ay*dt^2;
      z(k) = z(k-1) + vz(k-1)*dt + .5*az*dt^2;               
    end
    plot(y,z,'k');hold on
    plot(y(1),z(1),'go','MarkerFaceColor','green','markersize',12)
    plot(y(end),z(end),'ro','MarkerFaceColor','red','markersize',12)
    % labels for 2D plot       
    xlabel('y');ylabel('z');
    ttl=sprintf(['v_{z0}=' num2str(vz(1)) ' green = start, red = end, tstop = ' num2str(nt*dt) ', q=' num2str(q) ]);
    title(ttl)
    box
    hca = gca;
    xlim = hca.XLim;
    ylim = hca.YLim;
    [Y,Z] = meshgrid(linspace(xlim(1),xlim(2),50),linspace(ylim(1),ylim(2),50));
    B = Y;
    hold(hca,'on')
    %[c,h] = contour(Y,Z,B);    
    [h] = pcolor(Y,Z,B); 
    h.FaceAlpha = 0.3;
    shading flat;
    
    %clabel(c,h)
    colorbar
    cmap = irf_colormap('poynting');
    colormap(cmap)
    hold(hca,'off')
    hca.CLim = max(abs(hca.CLim))*[-1 1];
    
end
set(gcf,'position',[560   391   350   557])
axis equal
hca.XLim = xlim;
hca.YLim = ylim;
    
%% Speiser motions: Magnetosheath boundary layer



% Fields
tint = irf.tint('2015-10-16T10:33:30.20Z/2015-10-16T10:33:30.5Z'); 
time = gseE1.tlim(tint).time; % choose timeline
tN = double(time.ttns-time.start.ttns)*1e-9; % timeline in seconds

vN = -55*[-0.90 -0.28 -0.33]; % normal velocity, km/s, GSE
xN = tN*norm(vN); % normal length, km

BM = -gseB1.tlim(tint).resample(time).abs;
EN = gseE1.tlim(tint).resample(time).dot(vN/norm(vN));

% Fit of electric field
nPol = 6;
P = polyfit(xN,EN.data,nPol);
funEN = @(x,y,z) P(1)*x.^nPol + P(2)*x.^(nPol-1) + P(3)*x.^(nPol-2) + P(4)*x.^(nPol-3) + P(5)*x.^(nPol-4) + P(6)*x.^(nPol-5) + P(7)*x.^(nPol-6);
%funEN = @(x) -2 + 1e-3*x.^2.*abs(x-30);%x.*(100-2*x);%.*(x-10).^2;
% Fit of magnetic field
nPol = 6;
P = polyfit(xN,BM.data,nPol);
funBM = @(x,y,z) P(1)*x.^nPol + P(2)*x.^(nPol-1) + P(3)*x.^(nPol-2) + P(4)*x.^(nPol-3) + P(5)*x.^(nPol-4) + P(6)*x.^(nPol-5) + P(7)*x.^(nPol-6);

plot(xN,EN.data,xN,funEN(xN),xN,BM.data,xN,funBM(xN))
%%
close
% Constants
q=-units.e;m=units.me;
nt=29000*4;
dt=1e-11;

% initial conditions
x0=0;y0=5;z0=3;
r0=[x0,y0,z0]; % position vector

vz0s=[-5000 0.1 1 2 3 4 9 0.01 .5 1 2 4 10]*1e3;
vy0s=[0 0 0 0 0 0 0 0 0 0];
vx0s=[0 0 0 0 0 0 0 0 0 0];

% Set up fields for integration
% xyz = lmn
EX = @(x,y,z) 0;
EY = @(x,y,z) 0;
EZ = funEN;

BX = @(x,y,z) 0;
BY = funBM;
BZ = @(x,y,z) 0;

% intialise variables
x = zeros(nt,1); x(1) = x0;
y = zeros(nt,1); y(1) = y0;
z = zeros(nt,1); z(1) = z0;

for ii=[1]     
    vx = zeros(nt,1); vx(1) = vx0s(ii);
    vy = zeros(nt,1); vy(1) = vy0s(ii);
    vz = zeros(nt,1); vz(1) = vz0s(ii);
    
    for k=2:(nt)
      Ex = EX(x(k-1),y(k-1),z(k-1))*1e-3;
      Ey = EY(x(k-1),y(k-1),z(k-1))*1e-3;
      Ez = EZ(x(k-1),y(k-1),z(k-1))*1e-3;
      
      Bx = BX(x(k-1),y(k-1),z(k-1))*1e-9;
      By = BY(x(k-1),y(k-1),z(k-1))*1e-9;
      Bz = BZ(x(k-1),y(k-1),z(k-1))*1e-9;
      
      
      % Acceleration, depend on velocity, magnetic and electric field
      ax = q/m*(Ex+vy(k-1)*Bz-vz(k-1)*By);
      ay = q/m*(Ey+vz(k-1)*Bx-vx(k-1)*Bz);
      az = q/m*(Ez+vx(k-1)*By-vy(k-1)*Bx);
      
      % Velocities
      vx(k) = vx(k-1) + ax*dt;
      vy(k) = vy(k-1) + ay*dt;
      vz(k) = vz(k-1) + az*dt;

      % Positions
      x(k) = x(k-1) + vx(k-1)*dt + .5*ax*dt^2;
      y(k) = y(k-1) + vy(k-1)*dt + .5*ay*dt^2;
      z(k) = z(k-1) + vz(k-1)*dt + .5*az*dt^2;               
    end
    if 0
      plot3(x,y,z,'k');hold on
      plot3(x(1),y(1),z(1),'go','MarkerFaceColor','green','markersize',12)
      plot3(x(end),y(end),z(end),'ro','MarkerFaceColor','red','markersize',12)            
      % labels for 2D plot       
      xlabel('y');ylabel('z');
    elseif 0
      plot3(z,-y,x,'k');hold on
      plot3(z(1),-y(1),x(1),'go','MarkerFaceColor','green','markersize',12)
      plot3(z(end),-y(end),x(end),'ro','MarkerFaceColor','red','markersize',12)            
      % labels for 2D plot       
      xlabel('N');ylabel('M');zlabel('L');
      hca = gca;
      hca.XDir = 'reverse';
    else
      plot(z,x,'k');hold on
      plot(z(1),x(1),'go','MarkerFaceColor','green','markersize',12)
      plot(z(end),x(end),'ro','MarkerFaceColor','red','markersize',12)            
      % labels for 2D plot       
      xlabel('N');zlabel('L');
      hca = gca;
      hca.XDir = 'reverse';
      if 1 
        xlim = hca.XLim;
        ylim = hca.YLim;
        [X,Y] = meshgrid(linspace(xN(1),xN(end),20),linspace(ylim(1),ylim(2),20));
        [X,Y] = meshgrid(linspace(xlim(1),xlim(end),20),linspace(ylim(1),ylim(2),20));
        B = BY(X,Y,Y*0);
        hold(hca,'on')
        %[c,h] = contour(Y,Z,B);    
        [h] = contour(X,Y,B); 
        h.FaceAlpha = 0.3;
        shading flat;

        %clabel(c,h)
        colorbar
        %cmap = irf_colormap('poynting');
        %colormap(cmap)
        hold(hca,'off')
        %hca.CLim = max(abs(hca.CLim))*[-1 1];
      end
    end
    ttl=sprintf(['v_{z0}=' num2str(vz(1)) ' green = start, red = end, tstop = ' num2str(nt*dt) ', q=' num2str(q) ]);
    title(ttl)
    box
    
    if 0
      hca = gca;
      xlim = hca.XLim;
      ylim = hca.YLim;

      [Y,Z] = meshgrid(linspace(xlim(1),xlim(2),50),linspace(ylim(1),ylim(2),50));
      B = Y;
      hold(hca,'on')
      %[c,h] = contour(Y,Z,B);    
      [h] = pcolor(Y,Z,B); 
      h.FaceAlpha = 0.3;
      shading flat;

      %clabel(c,h)
      colorbar
      cmap = irf_colormap('poynting');
      colormap(cmap)
      hold(hca,'off')
      hca.CLim = max(abs(hca.CLim))*[-1 1];
    end
end
%%
set(gcf,'position',[560   391   350   557])
axis equal
hca.XLim = xlim;
hca.YLim = ylim;





    