b=-1:1;
b=[b;b;b];
a=permute(b,[2 1])-1;
c=[1 0; 0 -1];
h=irf_plot(1);
surfc(h,a,b,a*0,c)
hold(h,'on')
xlabel('x')
ylabel('y')
set(h,'ylim',[-2 2],'xlim',[-2 2])
%%
surf(h,-flipdim(a,2),b,a*0,c)