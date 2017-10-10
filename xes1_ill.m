% xes1


z = @(k,x,o,t,v0) cos(k*x-o*t-v0*k*t) + cos(k*x+o*t-v0*k*t);

L = 10;
x = linspace(0,L,1000);
T = 1;
o = 2*pi/T;
k = 2*pi/(L/5);
nt =200; t = linspace(0,3*T,nt);
v0=o/k;
for ii=1:5:nt; 
    plot(x,z(k,x,o,t(ii),0))
    set(gca,'ylim',2*[-1 1])
    pause(0.1)
end

%%
x=linspace(0,5,100);
funct=@(x)(-sqrt(pi)*x.^3.*exp(-x.^2));

%%

gamma = @(ope,or,k,ve) (-sqrt(pi/8)*ope*or.^2./k^3./ve^3).*exp(-(or/k/ve).^2/2);

ope = 1;
k   = 1.7;
ve  = 0.5;
%ve = 2*ve;;
nor=100;

or1=1.6

or = linspace(0,2,nor);
plot(or,gamma(ope,or,k,ve),...
        or1,gamma(ope,or1,k,ve),'ro',...
        or,gamma(ope,ones(nor)*or1,k,ve),'-r');
xlabel('\omega_r/\omega_p','fontsize',16)
ylabel('\gamma/\omega_p','fontsize',16)
ylim = get(gca,'ylim');
set(gca,'ylim',[ylim(1) 0.1])

%%

gamma = @(o,x) (-sqrt(pi)*x^3*o.^2).*exp(-(o*x).^2);

ope = 1;
k   = 1;
ve  = 0.5;
x=ope/k/ve;
nor=100;
o = linspace(0,2,nor);


plot(or,gamma(o,x)),...
        %1.3,gamma(o,x),'ro',...
        %or,gamma(ope,ones(nor)*1.3,k,ve),'-r',...
        %1,gamma(ope,1,k,ve),'go',...
        %or,gamma(ope,ones(nor)*1,k,ve),'-g');
xlabel('\omega_r/\omega_p','fontsize',16)
ylabel('\gamma/\omega_p','fontsize',16)
set(gca,'ylim',[-1.4 0.1])