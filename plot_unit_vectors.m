figure;
r=[1 0 0];
b=[0 0 1];
k=[0 0 0];
% normal along y
plot3([0 -0.83],[0 0.27],[0 0.49],'color',r,'linestyle',':'); hold(gca,'on')
plot3([0 -0.11],[0 -0.94],[0 0.33],'color',b,'linestyle',':'); hold(gca,'on')
plot3([0 0.55],[0 0.22],[0 0.81],'color',k,'linestyle',':'); hold(gca,'on')
% normal along x
plot3([0 0.16],[0 -0.97],[0 0.15],'color',b); hold(gca,'on')
plot3([0 0.82],[0 0.04],[0 -0.57],'color',r); hold(gca,'on')
plot3([0 0.55],[0 0.22],[0 0.81],'color',k); hold(gca,'on')

legend('normal 1','prop 1','z','prop 1','normal 2','')
axis equal
%%
figure;
r=[1 0 0];
b=[0 0 1];
k=[0 0 0];
% tool
plot3([0 -0.78],[0 0.62],[0 0.09],'color',b,'linestyle',':'); hold(gca,'on')
plot3([0 -0.28],[0 -0.48],[0 0.83],'color',r,'linestyle',':'); hold(gca,'on')
plot3([0 0.55],[0 0.62],[0 0.55],'color',k,'linestyle',':'); hold(gca,'on')
% mva
plot3([0 -0.82],[0 0.50],[0 0.27],'color',b); hold(gca,'on')
plot3([0 -0.11],[0 -0.60],[0 0.79],'color',r); hold(gca,'on')
plot3([0 0.55],[0 0.62],[0 0.55],'color',k); hold(gca,'on')

legend('prop 1','normal 1','z','prop 1','normal 2','')
axis equal