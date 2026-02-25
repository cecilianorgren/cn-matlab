[A,im,map]=irf_plasma_wave_visualization('demo','alfven_compr');
movie(A,20);
%%
E = [0 0.01 0];
kx = 1/0.5;
kz = 0;
f=0.4/(2*pi);
T=1/f;
N=400;
[A,im,map]=irf_plasma_wave_visualization('E',E,'kx',kx,...
    'kz',kz,'view','3d','N',N,'f',f,...
    't',1/f);
imwrite(im,map,'AC.gif','DelayTime',0,'LoopCount',inf);
movie(A,10);
%%
[A,im,map]=irf_plasma_wave_visualization('demo','lowerhybrid');
movie(A,20);
%%
lx=0.2;kx=20;%2*pi/lx;
lz=2;kz=4;%2*pi/lz;
E=[2 0 2*kz/kx];
qm=[1 -400];
f=3;%20/(2*pi);
N=500;
[A,im,map]=irf_plasma_wave_visualization('E',E,'kx',kx,...
    'kz',kz,'qm',qm,'view','xz','N',500,'f',f,...
    't',1/f,'z',1,'x',1);
imwrite(im,map,'LH.gif','DelayTime',0,'LoopCount',inf);
movie(A,10);


%%
[A,im,map]=irf_plasma_wave_visualization('E',[0 1 0.01],'ky',2*2*pi,...
    'kx',0,'kz',2*0.2*pi,'N',500,'f',6,'T',10);
movie(A,10);
