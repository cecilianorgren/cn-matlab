%% FPI, 1 DES

nPolar = 16;  polar = linspace(-90,90,nPolar+1);
nAzimuth = 4; azim  = linspace(5,360/4-5,nAzimuth+1);

[POL,AZ] = meshgrid(polar,azim);
X = cosd(POL).*cosd(AZ);
Y = cosd(POL).*sind(AZ);
Z = sind(POL);
C = AZ;
%for ic = 1:numel(C); C(ic) = 1; end
surf(X,Y,Z,C)
cmap = pic_colors('thermal');
colormap(cmap(20:end,:))
cmap = [1 0 0; 0 1 0; 0 1 1; 1 0 1]*0.95;
colormap(cmap)
axis equal;
axis off
%% FPI, 4 DES
nPolar = 16;  polar = linspace(-90,90,nPolar+1);
nAzimuth = 8; azim  = linspace(2,360/nAzimuth-2,nAzimuth+1);

plot(nan,nan)
hold on
for ides = 1:8
  az_shift = 45;
  azim_tmp = azim + (ides-1)*az_shift;
  [POL,AZ] = meshgrid(polar,azim_tmp);
  X = cosd(POL).*cosd(AZ);
  Y = cosd(POL).*sind(AZ);
  Z = sind(POL);
  C = AZ - (ides-1)*az_shift;
  %for ic = 1:numel(C); C(ic) = 1; end
  surf(X,Y,Z,C)
  
  cmap = pic_colors('thermal');
  colormap(cmap(20:end,:))
  axis equal;
end
hold off
axis off

%% FPI, 1 DES, 2 energies

nPolar = 16;  polar = linspace(-90,90,nPolar+1);
nAzimuth = 4; azim  = linspace(5,360/4-5,nAzimuth+1);

[POL,AZ] = meshgrid(polar,azim);
X = cosd(POL).*cosd(AZ);
Y = cosd(POL).*sind(AZ);
Z = sind(POL);
C = AZ;
%for ic = 1:numel(C); C(ic) = 1; end
surf(X,Y,Z,C)
cmap = pic_colors('thermal');
colormap(cmap(20:end,:))
cmap = [1 0 0; 0 1 0; 0 1 1; 1 0 1]*0.95;
colormap(cmap)
axis equal;

hold on
R = 1.25;
[POL,AZ] = meshgrid(polar,azim);
X = R*cosd(POL).*cosd(AZ);
Y = R*cosd(POL).*sind(AZ);
Z = R*sind(POL);
C = AZ+0;
hs = surf(X,Y,Z,C);
hs.FaceAlpha = 1;
hold off


hold on
R = 1.5;
[POL,AZ] = meshgrid(polar,azim);
X = R*cosd(POL).*cosd(AZ);
Y = R*cosd(POL).*sind(AZ);
Z = R*sind(POL);
C = AZ+0;
hs = surf(X,Y,Z,C);
hs.FaceAlpha = 1;
hold off

axis off