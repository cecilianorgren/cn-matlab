% Plot illustration of deHoffmanTeller frame
vA = 4;
Bn = 1;
Bt = 7;
B1 = [Bt  Bn]; % [Bt Bn]
B2 = [-Bt Bn];
v1_HT = vA*B1/norm(B1);
v2_HT = vA*B2/norm(B2);
vHT = [2 -1];
v1 = v1_HT + vHT;
v2 = v2_HT + vHT;

colors = mms_colors('matlab');
fontsize = 15;

hca = subplot(1,1,1);
plot(hca,[-Bt 1]*1.2,[0 0],'k--'); 
hold(hca,'on')

quiver(hca,-8,0,0,1,'k')
text(-7.7,0.6,'n','color',[0 0 0],'fontsize',fontsize)

quiver(hca,-B1(1),-B1(2),B1(1),B1(2),1,'k')
quiver(hca,0,0,B2(1),B2(2),1,'k')
text(-6.2,0.6,'B','color',[0 0 0],'fontsize',fontsize)

quiver(hca,-v1_HT(1),-v1_HT(2),v1_HT(1),v1_HT(2),1,'color',colors(5,:))
quiver(hca,0,0,v2_HT(1),v2_HT(2),1,'color',colors(5,:))
text(-3,-0.2,'v''_{1}','color',colors(5,:),'fontsize',fontsize)
text(-2.5,0.5,'v''_{2}','color',colors(5,:),'fontsize',fontsize)

quiver(hca,-v1(1),-v1(2),v1(1),v1(2),1,'color',colors(2,:))
quiver(hca,0,0,v2(1),v2(2),1,'color',colors(2,:))
text(-2,-0.7,'v_{1}','color',colors(2,:),'fontsize',fontsize)
text(-1.2,1.2,'v_{2}','color',colors(2,:),'fontsize',fontsize)

quiver(hca,-v1(1),-v1(2),vHT(1),vHT(2),1,'color',colors(3,:))
quiver(hca,v2_HT(1),v2_HT(2),vHT(1),vHT(2),1,'color',colors(3,:))
text(-6,-1.1,'v_{HT}','color',colors(3,:),'fontsize',fontsize)
text(-4,1.1,'v_{HT}','color',colors(3,:),'fontsize',fontsize)

quiver(hca,-v1(1)+vHT(1)  ,-v1(2)  ,0,vHT(2),1,'color',colors(3,:),'linestyle','--')
quiver(hca,v2_HT(1)+vHT(1),v2_HT(2),0,vHT(2),1,'color',colors(3,:),'linestyle','--')
text(-6,-1.1,'v_{HT}','color',colors(3,:),'fontsize',fontsize)
text(-4,1.1,'v_{HT}','color',colors(3,:),'fontsize',fontsize)

text(-3.8,-1.3,'u_n = n\cdot v_{HT}','color',colors(3,:),'fontsize',fontsize)
text(-1.7,0.8,'u_n = n\cdot v_{HT}','color',colors(3,:),'fontsize',fontsize)


hold(hca,'off')

box(hca,'off')
axis(hca,'off')


