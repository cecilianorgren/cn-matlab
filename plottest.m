%% Combining plots
figure(50)
plot(0:.01:1)
figure(60)
plot((0:.01:1).^2)


set(gca,'YDir','reverse');

L = findobj(50,'type','line');
copyobj(L,findobj(60,'type','axes'));