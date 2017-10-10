% Make semitransparent positionplot
global rlim
nParticles = 50000;
%
r0 = rlim*randn(nParticles,1); 
%r0(find(abs(r0)>rlim)) = [];
nParticles= numel(r0)
angle = 360*rand(nParticles,1);
x0 = r0.*cosd(angle);
y0 = r0.*sind(angle);


m = log((r0/rlim).^2-1);%./r0;
m=ones(nParticles,1);
m = abs(r0).*rlim*2*pi.*exp(r0.^2/rlim.^2/2);
%m = m/max(m);

%%
r0 = rlim*rand(nParticles,1);
angle = 360*rand(nParticles,1);
x0 = r0.*cosd(angle);
y0 = r0.*sind(angle);

m = r0/rlim;%./r0;
%m = 

%%
sh = scatter(x0,y0,m*100);
h = findobj(sh,'Type','patch');
set(h,'facealpha',0.5);
%%
figure(2)
circle = 0:pi/10:2*pi;
for ip=1:nParticles
    hs=patch(sin(circle)+x0(ip),cos(circle)+y0(ip),'b','edgecolor','none');
    alpha(hs,m(ip)*0.1)
end

%%
circle = 0:pi/10:2*pi;
for ip=1:nParticles
    hs=patch(sin(t)+x0(ip)*m0(ip),cos(t)+y0(ip)*m0(ip),'b','edgecolor','none');
    alpha(hs,0.1)
end

%%
circle = 0:pi/10:2*pi;
for ip=1:nParticles
    hs=patch(sin(t)+m0(ip),cos(t)+y0(ip)*m0(ip),'b','edgecolor','none');
    alpha(hs,0.1)
end

%%
circle = 0:pi/10:2*pi;
for ip=1:nParticles
    hs=patch(sin(t)+x0(ip),cos(t)+y0(ip),'b','edgecolor','none');
    alpha(hs,0.1)
end
%% Example
x=randn(5000,1)*20;
y= randn(5000,1)*20;
t= 0:pi/10:2*pi;
figure();
for i=1:size(x0)
    pb=patch((sin(t)+ x0(i)),(cos(t)+y0(i)),'b','edgecolor','none');
    alpha(pb,.1);
end
