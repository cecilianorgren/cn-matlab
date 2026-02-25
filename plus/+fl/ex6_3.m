%ex6.3
v0=3;
R1=1;
R2=3;
v = @(r)(v0*log(r/R1)/log(R2/R1));

r=linspace(R1,R2,40);

plot(v(r)/v0,r)
ylabel('r')
xlabel('v(r)/v_0')

plot(r,v(r)/v0)
xlabel('r')
ylabel('v/v_0')