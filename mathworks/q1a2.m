clc
clf
clear all

% % Author
%  Judah S
% Sixth Order RK-Fehlberg Equation
% Function y' + 2xy = 0; initial condition y(0)=1; where 0 <= x <= 10
% You can easily change the function and rest condition according to your
% requirment.
   fcnstr='-2*x*y' ;
   f=inline(fcnstr) ;

   x0=0 ;
   y0=1 ;

% xf, x location at where you wish to see the solution to the ODE

   xf=10 ;

% n, number of steps to take

   n=5 ;

format long g
h=(xf-x0)/n ;
xa(1)=x0 ;
ya(1)=y0 ;

for i=1:n
    % Adding Step Size
  xa(i+1)=xa(i)+h ;

% Calculating k1, k2, k3, k4, K5 and K6
  k1 = f(xa(i),ya(i)) ;
  k2 = f(xa(i)+(1/4)*h,ya(i)+(1/4)*k1*h) ;
  k3 = f(xa(i)+(3/8)*h,ya(i)+(3/32)*h*(k1+3*k2)) ;
  k4 = f(xa(i)+(12/13)*h,ya(i)+(12/2197)*h*(161*k1-600*k2+608*k3)) ;
  k5 = f(xa(i)+h,ya(i)+(1/4104)*h*(8341*k1-32832*k2+29440*k3-845*k4)) ;
  k6 = f(xa(i)+(0.5)*h,ya(i)+h*(-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5)) ;
% Using 6th Order Runge-Kutta formula
  ya(i+1)=ya(i)+1/5*((16/27)*k1+(6656/2565)*k3+(28561/11286)*k4-(9/10)*k5+(2/11)*k6)*h ;

  end


disp(sprintf('\n\n********************************Results**********************************'))

% The following finds what is called the 'Exact' solution
xspan = [x0 xf];
[x,y]=ode45(f,xspan,y0);
[yfi dummy]=size(y);
yf=y(yfi);

% Plotting the Exact and Approximate solution of the ODE.
hold on
xlabel('x');ylabel('y');
title('Exact and Approximate Solution of the ODE by the 4th Order Runge-Kutta Method');
plot(x,y,'--','LineWidth',2,'Color',[0 0 1]);            
plot(xa,ya,'-','LineWidth',2,'Color',[0 1 0]);
legend('Exact','Approximation');


disp(sprintf('\n   Approximate = %g',ya(n+1))) 
disp(sprintf('   Exact       = %g',yf))
disp(sprintf('\n   True Error = Exact - Approximate')) 
disp(sprintf('              = %g - %g',yf,ya(n+1)))
disp(sprintf('              = %g',yf-ya(n+1)))
disp(sprintf('\n   Absolute Relative True Error Percentage'))
disp(sprintf('              = | ( Exact - Approximate ) / Exact | * 100'))
disp(sprintf('              = | %g / %g | * 100',yf-ya(n+1),yf))
disp(sprintf('              = %g',abs( (yf-ya(n+1))/yf )*100))

