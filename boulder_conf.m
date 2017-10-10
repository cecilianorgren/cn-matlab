figure(74)
setupfigure
set(gcf,'defaultAxesFontSize',16);
set(gcf,'defaultTextFontSize',19);
plot(0,0,'bo','markersize',15,'linewidth',3); 
hold(gca,'on'); 
plot(-facdy,-facdz,'go','markersize',15,'linewidth',3);
text(0,0,'   C4');
text(-facdy,-facdz,'   C3');
set(gca,'ylim',[-12 48],'xlim',[-8 23])
axis square;
xlabel('y - measured component in sc spin plane')
ylabel('z - estimated component parallel to B')
arrow([15 0],[0 20]+[15 0]); text(17,3,'B');
set(gcf,'PaperPositionMode','auto');
%%
s=0.5;
%quiver(0,0,-10*s,0)
arrow([0 0],[-12 0]*s);arrow([0 0],[0 -19]*s);arrow([0 0],[0 20]*s)
arrow(-[facdy facdz],[-20 0]*s-[facdy facdz]);
arrow(-[facdy facdz],[0 -20]*s-[facdy facdz]); 
arrow(-[facdy facdz],[0 20]*s-[facdy facdz]);