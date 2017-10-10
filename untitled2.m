bbB4=[bB4(:,1) bB4(:,2:4)*0.4];
bbB3=[bB3(:,1) bB3(:,2:4)*0.4];
figure;irf_plot(phi3(:,1:2),'color',[0.8 0.8 0.2]);
hold on;
irf_plot(bbB3(:,[1 4]),'r');