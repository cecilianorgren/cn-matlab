figure;x=(0:0.01:1)*i;
P=@(x,rat)(1-(rat./x).^2);
kvr=@(x,rat)(x.^2+1./P(x,rat));
kvi=@(x,rat)(sqrt(-(1+4*P(x,rat).*x.^2))./P(x,rat));
ratio=0.6;
plot(kvr(x,ratio),kvi(x,ratio),kvr(x,ratio),-kvi(x,ratio));hold on;