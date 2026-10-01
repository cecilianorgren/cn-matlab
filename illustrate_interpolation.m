t = 1;
ffun = @(x,t) exp(-(x).^2/t^2);

%x_edges = linspace(0,1,11);
x_edges = logspace(-2,log10(5),16);
x_centers = x_edges(1:end-1) + 0.5*diff(x_edges);

f1 = ffun(x_centers-0.0,t);
f2 = ffun(x_centers-2.0,t);
xq = logspace(-2,2,72);

method = 'linear';
fint_1 = interp1(x_centers,log10(f1),xq,method);
fint_1 = 10.^fint_1;

fint_2 = interp1(x_centers,log10(f2),xq,method);
fint_2 = 10.^fint_2;

fint_3 = interp1(x_centers,(f1),xq,method);
%fint_1 = 10.^fint_1;

fint_4 = interp1(x_centers,(f2),xq,method);
%fint_2 = 10.^fint_2;


hca = subplot(4,1,1);
loglog(hca,x_centers,f1,'-*',xq,fint_1,'o')

hca = subplot(4,1,2);
loglog(hca,x_centers,f2,'-*',xq,fint_2,'o')

hca = subplot(4,1,3);
loglog(hca,x_centers,f1,'-*',xq,fint_3,'o')

hca = subplot(4,1,4);
loglog(hca,x_centers,f2,'-*',xq,fint_4,'o')
