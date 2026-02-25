% fluid ex 6.8


y = 0:0.05:1;
vx1 = @(y)(y-y.^2/2);
plot(vx1(y),y)

title('Exercise 6.8 - normalized velocity')
xlabel('v_x')
ylabel('y')

%% fluid ex 6.9

y = 0:0.05:1;
vx2 = @(y)(y-y.^2);
plot(vx2(y),y,vx1(y),y)

title('Exercise 6.9 & 6.9 - normalized velocity')
xlabel('v_x')
ylabel('y')

%% fluid ex 6.8 & 6.9

y = 0:0.05:1;
vx1 = @(y) y-y.^2/2;
vx2 = @(y) (y-y.^2)/2;
lines =plot(vx2(y),y,vx1(y),y);

title('Normalized velocity profiles')
xlabel('v_x/(\rho g sin{\theta}/\eta)')
ylabel('y/h')
%legend('solid surface at y=h','free surface at y=h','location','southeast')
text(0.12,0.8,{'solid fix surface at y=h','v_x = 0.5 (yh - y^2)'},'fontsize',15,'color',lines(1).Color)
text(0.3,0.25,{'free surface at y=h','v_x = yh - 0.5y^2'},'fontsize',15,'color',lines(2).Color)
hca = gca;
hca.FontSize = 15;