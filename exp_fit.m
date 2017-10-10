p1=polyfit(specrec.en,log(PA0),1);
p2=polyfit(specrec.en,log(PA0),2);

f1=@(x)(exp(p1(2)).*exp(p1(1)*x));
f2=@(x)(exp(p2(3)).*exp(p2(2)*x.^1).*exp(p2(1)*x.^2));

loglog(specrec.en,PA0,specrec.en,f1(specrec.en),specrec.en,f2(specrec.en))
legend('Data','e^de^{fx}','e^ae^{bx}he^{cx^2}')
%%
p1=polyfit(specrec.en,log(PA90),1);
p2=polyfit(specrec.en,log(PA90),2);

f1=@(x)(exp(p1(2)).*exp(p1(1)*x));
f2=@(x)(exp(p2(3)).*exp(p2(2)*x.^1).*exp(p2(1)*x.^2));

loglog(specrec.en,PA90,specrec.en,f1(specrec.en),specrec.en,f2(specrec.en))
legend('Data','e^de^{fx}','e^ae^{bx}he^{cx^2}')
%%
i0=1:length(x);
p1=polyfit(x(i0),log(y2(i0)),2);
f1=@(x)(exp(p1(3)).*exp(p1(2)*x).*exp(p1(1)*x.^2));

i1=8:length(x);
p2=polyfit(x(i1),log(y2(i1)),2);
f2=@(x)(exp(p2(3)).*exp(p2(2)*x.^1).*exp(p2(1)*x.^2));

loglog(x,y2,x,f1(x),x,f2(x),x,f1(x)+f2(x),x(i0),y2(i0),'o')
legend('Data','f0','f1','f0+f1','f0 data')
%%
p0=polyfit(x,log(y2),2);
f0=@(x)(exp(p0(3)).*exp(p0(1)*x).*exp(p0(1)*x.^2));

yrest=y2-f0(x);
p1=polyfit(x(8:12),log(yrest(8:12)),2);
f1=@(x)(exp(p1(3)).*exp(p1(1)*x).*exp(p1(1)*x.^2));

plot(x,y2,x,f0(x),x,f1(x))
legend('Data','e^de^{fx}','e^ae^{bx}he^{cx^2}')
%%
p1=polyfit(specrec.en,log(PA180),1);
p2=polyfit(specrec.en,log(PA180),2);

f1=@(x)(exp(p1(2)).*exp(p1(1)*x));
f2=@(x)(exp(p2(3)).*exp(p2(2)*x.^1).*exp(p2(1)*x.^2));

plot(specrec.en,PA180,specrec.en,f1(specrec.en),specrec.en,f2(specrec.en))
legend('Data','e^de^{fx}','e^ae^{bx}he^{cx^2}')
%%
r=-3:0.01:3;
t0=1*exp(-(r).^2/10);
t1=exp(-(r-1).^2/0.05);
plot(r,t1)
%%
r=-3:0.01:3;
t0=10*exp(-(x).^2/(1e2)^2);
t1=exp(-(x-0.1e4).^2/(1e2)^2);
plot(x,t0,x,t1,x,t0+t0)
legend('t0','t1','t0+t1')
%% Polyfit manuellt s? man kan ha std ocks?
G=[sum(x.^2) sum(x); sum(x) length(x)];
invG=G^-1;
h=[sum(x.*log(y));sum(log(y))];
a=invG*h;
f=exp(a(1)*x)*exp(a(2));
plot(x,y,'*',x,f)
S=sum((f-log(y)).^2)