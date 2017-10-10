nx = 1000;
L = 1;
x = linspace(0,2*L,nx);
om = 1; T = 2*pi/om;
k = 2*pi/L;
nt = 100;
t = linspace(0,T,nt);
dt = t(2)-t(1);
g = 10; 

% shallow waves
H = L/10;
om = @(k,g) sqrt(g*H)*k;
vp_shal = @(k,g) sqrt(g*H);

% deep waves
H = 10*L;
om = @(k,g) sqrt(g*k);
vp_deep = @(k,g,o) g./o;

% the wave
v0 = -vp_shal(k,g);
v0 = -vp_deep(k,g,om(k,g));
wave = @(a,k,x,om,t,v0) a*sin(k*x-om*t-v0*k*t);

%v0=0;
for ii = 1:nt
    plot(x,wave(1,k,x,om(k,g),t(ii),v0));
    
    %totwave = wave(a(1),k(1),x,om(k(1),g),t(ii),v0) + wave(a(2),k(2),x,om(k(2),g),t(ii),v0) + wave(a(3),k(3),x,om(k(3),g),t(ii),v0);
    %subplot(4,1,4); plot(x,totwave); hold on;
    %plot(0.5*l(1)*[1 1],sum(a)*[-1 1],col{1})
    %plot(0.5*l(2)*[1 1],sum(a)*[-1 1],col{2})
    %plot(0.5*l(3)*[1 1],sum(a)*[-1 1],col{3})
    %plot(x,x*0)
    %plot(vg(g,k(1))*t(ii),0,'bo');plot(vp(g,k(1))*t(ii),0,'bx');
    %plot(vg(g,k(2))*t(ii),0,'go');plot(vp(g,k(2))*t(ii),0,'gx');
    %plot(vg(g,k(3))*t(ii),0,'ro');plot(vp(g,k(3))*t(ii),0,'rx');
    %set(gca,'ylim',sum(a)*[-1 1],'xlim',x([1 end]))
    hold off;
    title(['t = ' num2str(t(ii))])
    drawnow
    pause(0.1)
end


%%


nx = 1000; 
x = linspace(0,2*2*pi,nx);

l = 2*pi*[0.95 1 1.05];
k = 2*pi./l;
k = [1.05 1 0.95];
k1 = 1.05;
k2 = 1;
k3 = 2*k1*k2/(k1+k2);

k=[k1 k2 k3];
k=[0.99 1 1.01];
nt = 50;
t = linspace(0,2000,nt);

g = 10; 
om = @(k,g) sqrt(g*k);

vg = @(g,k) sqrt(g./k)/2;
vp = @(g,k) sqrt(g./k);

wave = @(a,k,x,om,t,v0) a*sin(k*x-om*t-v0*k*t);
a = [0,1,1];
include = find(a~=0);
col = {'b','g','r'};
%k0=2*pi/mean(l(include));
%k0 =  median(k);
%v0 = -vp(g,k0);
v0 = -vp(g,mean(k(include)));
for ii = 1:nt
    for jj = 1:numel(k)
        subplot(4,1,jj); plot(x,wave(a(jj),k(jj),x,om(k(jj),g),t(ii),v0),col{jj});
    end
    totwave = wave(a(1),k(1),x,om(k(1),g),t(ii),v0) + wave(a(2),k(2),x,om(k(2),g),t(ii),v0) + wave(a(3),k(3),x,om(k(3),g),t(ii),v0);
    subplot(4,1,4); plot(x,totwave); hold on;
    plot(0.5*l(1)*[1 1],sum(a)*[-1 1],col{1})
    plot(0.5*l(2)*[1 1],sum(a)*[-1 1],col{2})
    plot(0.5*l(3)*[1 1],sum(a)*[-1 1],col{3})
    plot(x,x*0)
    %plot(vg(g,k(1))*t(ii),0,'bo');plot(vp(g,k(1))*t(ii),0,'bx');
    %plot(vg(g,k(2))*t(ii),0,'go');plot(vp(g,k(2))*t(ii),0,'gx');
    %plot(vg(g,k(3))*t(ii),0,'ro');plot(vp(g,k(3))*t(ii),0,'rx');
    set(gca,'ylim',sum(a)*[-1 1])
    hold off;
    pause(0.1)
end


%%
nx = 1000;
x = linspace(0,40*2*pi,nx);

l = 2*pi*[0.95 1 1.05];
k = 2*pi./l;
k = [1.05 1 0.95];
k1 = 1.05;
k2 = 1;
k3 = 2*k1*k2/(k1+k2);

k=[k1 k2 k3];
nt = 50;
t = linspace(0,600,nt);

g = 10; 

om = @(k,g,H) sqrt(g*H)*k;
om = @(k,g) om(k,g,H);
vg = @(g,H) sqrt(g*H);
vp = @(g,H) sqrt(g*H);


wave = @(a,k,x,om,t,v0) a*sin(k*x-om*t-v0*k*t);
a = [1,1,1];
include = find(a~=0);
col = {'b','g','r'};
k0=2*pi/mean(l(include));
k0 =  median(k);
v0 = -vp(g,H);
for ii = 1:nt
    for jj = 1:numel(k)
        subplot(4,1,jj); plot(x,wave(a(jj),k(jj),x,om(k(jj),g),t(ii),v0),col{jj});
    end
    totwave = wave(a(1),k(1),x,om(k(1),g),t(ii),v0) + wave(a(2),k(2),x,om(k(2),g),t(ii),v0) + wave(a(3),k(3),x,om(k(3),g),t(ii),v0);
    subplot(4,1,4); plot(x,totwave); hold on;
    plot(0.5*l(1)*[1 1],sum(a)*[-1 1],col{1})
    plot(0.5*l(2)*[1 1],sum(a)*[-1 1],col{2})
    plot(0.5*l(3)*[1 1],sum(a)*[-1 1],col{3})
    plot(x,x*0)
    %plot(vg(g,k(1))*t(ii),0,'bo');plot(vp(g,k(1))*t(ii),0,'bx');
    %plot(vg(g,k(2))*t(ii),0,'go');plot(vp(g,k(2))*t(ii),0,'gx');
    %plot(vg(g,k(3))*t(ii),0,'ro');plot(vp(g,k(3))*t(ii),0,'rx');
    set(gca,'ylim',sum(a)*[-1 1])
    hold off;
    pause(0.1)
end
%%

nx = 2000;
L = 1;
x = linspace(0,4*L,nx);
om = 1; T = 2*pi/om;
k = 2*pi/L;
nt = 100;
t = linspace(0,T/2,nt);
dt = t(2)-t(1);
g = 10; 
wave_type = 'deep';

switch wave_type
    case 'shallow' % shallow waves
        H = L/10;
        om = @(k,g) sqrt(g*H)*k;
        vp_shal = @(k,g) sqrt(g*H);
        v0 = -vp_shal(k,g);
    case 'deep' % deep waves
        H = 10*L;
        om = @(k,g) sqrt(g*k);
        vp_deep = @(k,g,o) g./o;
        v0 = -vp_deep(k,g,om(k,g));
end
% the wave
wave = @(k,x,om,t,v0) exp(1i*(k*x-k*v0*t-om*t)) + exp(-1i*(k*x-k*v0*t+om*t));

%v0=0;
v0 = -v0;
x0 = 0;
for ii = 1:nt
    x0 = x0+dt*v0;
    plot(x,real(wave(k,x,om(k,g),t(ii),v0)),x0*[1 1],2*[-1 1]);
    
    %totwave = wave(a(1),k(1),x,om(k(1),g),t(ii),v0) + wave(a(2),k(2),x,om(k(2),g),t(ii),v0) + wave(a(3),k(3),x,om(k(3),g),t(ii),v0);
    %subplot(4,1,4); plot(x,totwave); hold on;
    %plot(0.5*l(1)*[1 1],sum(a)*[-1 1],col{1})
    %plot(0.5*l(2)*[1 1],sum(a)*[-1 1],col{2})
    %plot(0.5*l(3)*[1 1],sum(a)*[-1 1],col{3})
    %plot(x,x*0)
    %plot(vg(g,k(1))*t(ii),0,'bo');plot(vp(g,k(1))*t(ii),0,'bx');
    %plot(vg(g,k(2))*t(ii),0,'go');plot(vp(g,k(2))*t(ii),0,'gx');
    %plot(vg(g,k(3))*t(ii),0,'ro');plot(vp(g,k(3))*t(ii),0,'rx');
    set(gca,'ylim',2*[-1 1],'xlim',x([1 end]))
    hold off;
    title(['t = ' num2str(t(ii))])
    drawnow
    pause(0.1)
end


