% extra exercise
eta = 1;
rho = 1;
V=-0:0.5:14;
alpha = V*rho/eta;
h = 1;
U=2;
y = 0:0.05:h;
[A,Y]=meshgrid(alpha,y);

u = @(alpha,y) U*(1-exp(alpha.*y))./(1-exp(alpha*h));
%U=
mesh(A,Y,u(A,Y))
xlabel('\alpha=V*\rho/\eta')
ylabel('y')
title('Velocity')

%%
y = linspace(0,h,100);
lines = plot([0 U]/U,[0 h],'k--',...
             u(0.1*rho*eta,y)/U,y,...
             u(2*rho*eta,y)/U,y,...
             u(3*rho*eta,y)/U,y,...
             u(4*rho*eta,y)/U,y,...
             u(10*rho*eta,y)/U,y,...
             u(40*rho*eta,y)/U,y);
xlabel('v_x/U')
%xlabel('\alpha=V*\rho/\eta')
ylabel('y/h')
title('Velocity profiles')
alphas = [0.1 2 3 4 10 40];
ys = [0.9 0.8 0.6 0.5 0.4 0.3];
ys = [0.42 0.55  0.6  0.65 0.8 0.9];
%ys = ones(1:10)*0.8;
for k = 1:numel(lines)-1
   %text(0.5,0.5,'s')
   %alphas(k)
   %u(alphas(k),0.5)
   %text(u(alphas(k),ys(k))/U,ys(k),[' \alpha = ' num2str(alphas(k))],'color',lines(k+1).Color,'fontsize',12)
   legends{k} = [' \alpha = ' num2str(alphas(k))];
end
hleg = legend({'\alpha \rightarrow 0',legends{:}},'location','southeast');
hleg.Box = 'Off';
hca=gca;
hca.FontSize = 15;
%text(0.78,0.55,'\alpha = \rhoV/\eta','fontsize',15)
text(0.82,0.55,'\alpha = V/\nu','fontsize',15)
if 0
hold on
quiver(0.5,0.1,0,0.08,3,'k')
text(0.5,0.2,'  V','fontsize',15)
hold off
end