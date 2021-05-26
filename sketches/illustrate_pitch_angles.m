r = [0 1];
theta = 0:11.25:360;
z = 0;
[R,TH] = meshgrid(r,theta);
[h,c] = polarPcolor(R,TH,R*0);

%%

R = [0 1];
theta = 0:11.25:180;
Z = (1:numel(theta))'*[1 1];
  figure
[h,c] = polarPcolor(R,theta,Z,'Ncircles',2,'Nspokes',numel(theta));