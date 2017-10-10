
x = linspace(0,pi,20);
y = x*0+1;
grad_n = cos(x)+1; 
plot3(x,y,grad_n,'linewidth',2,'color',[0 0 0]); hold on;

quiver3([0 0 0],[0 0 0],[0 0 0],[1 0 0],[0 1 0],[0 0 1],3,'linewidth',2,'color',[0 0 0])
quiver3([0.5 0.8],[0 3],[1 1],[0 0],1*[1 -1],[0 0],1.5,'linewidth',2,'color',[0 0 0])
%quiver3([0.8],[3],[1],[0],[-1],[0],1,'linewidth',2)
axis equal
view([1 1 0.7])
axis off
xlabel('x')
ylabel('y')
hold off

text()