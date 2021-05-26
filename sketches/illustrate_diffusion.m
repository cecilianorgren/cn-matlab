% illustrate diffusion

B = @(x,x0,l) tanh((x-x0)/l);
x = linspace(-5,5,1000);

%hca = subplot(1,2,1);
%plot(hca,x,B(x,0,1),'k-',x,B(x,1,1),'k--')

hca = subplot(1,1,1);
%plot(hca,x,B(x,0,1),'k-',x,B(x,0,2),'k--',x(2:end-1),diff(diff(B(x,0,1))),'r-')
plot(hca,x,B(x,0,1),'k-',x(2:end),diff(B(x,0,1))*50,'b-',x(2:end-1),diff(diff(B(x,0,1)))*2000,'r-',x(2:end-1),B(x(2:end-1),0,1)+diff(diff(B(x,0,1)))*2000,'k--')
legend(hca,{'B_1','\nabla B_1','\nabla^2B_1\propto \partial B_1/\partial t','B_2=B_1+\nabla^2B_1'},'location','southeast')