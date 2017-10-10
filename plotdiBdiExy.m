h=irf_plot(3)
irf_plot(h(1),diB4)
irf_plot(h(2),{diE3(end*3/4:end,[1 2]),diE4(end*3/4:end,[1 2])},'comp')
irf_plot(h(3),{diE3(end*3/4:end,[1 3]),diE4(end*3/4:end,[1 3])},'comp')
legend(h(1),'x','y','z','abs')
legend(h(2),'C3 x','C4 x')
legend(h(3),'C3 y','C4 y')

%%
h=irf_plot(3)
irf_plot(h(1),diB4)
irf_plot(h(2),{diE3(:,[1 2]),diE4(:,[1 2])},'comp')
irf_plot(h(3),{diE3(:,[1 3]),diE4(:,[1 3])},'comp')
legend(h(1),'x','y','z','abs')
legend(h(2),'C3 x','C4 x')
legend(h(3),'C3 y','C4 y')

