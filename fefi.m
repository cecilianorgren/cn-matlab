TiTe=3;
mime=1836;
di=0.01;
de=-di/TiTe;
x=-15:0.02:15;
veth=1;
vith=sqrt(TiTe/mime)*veth;

fe=(1-x*de/veth/veth).*exp(-(x-de).^2/(veth^2));
fi=(1-x*di/vith/vith).*exp(-(x-di).^2/(vith^2));

plot(x,fe,x,fi)
legend('e^-','i^+')
xlabel('v_y')
ylabel('f')
title('Ion and electron distribution functions')