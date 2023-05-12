bx = @(x,k,A,B) A*exp(k*x) + B*exp(-k*x);

a = 1;
x = linspace(-2,2,100);
k = 0.5;

C = 1;
A = C/(2*k*a)*exp(-2*k*a);
B = C/(2*k*a)*(2*k*a-1);
plot(x,bx(x,k,A,B))

%%

a = 1;
x = linspace(-2,2,100);
%k = 0.5;
k = 0.2:0.2:3.5;
[X,K] = ndgrid(x,k);
C = 1;

bx = tearing_bx(X,K,a,C);

plot(x,bx)

pcolor(x,k,bx')
colormap(pic_colors('blue_red'))
hb = colorbar(gca);