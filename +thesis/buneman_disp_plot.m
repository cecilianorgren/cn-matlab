% Buneman 3-comp
% run disprelESW_cn.m
% or
load /Users/Cecilia/MATLAB/cn-matlab/+thesis/dispcurv_modbun.mat
dispcurv1 = dispcurv;
load /Users/Cecilia/MATLAB/cn-matlab/+thesis/dispcurv_ee_s07.mat
dispcurv2 = dispcurv;

%
% Variables are defined as
%
% vith = sqrt(2*qe*Ti/mi);
% veth1 = sqrt(2*qe*Te1/me);
% veth2 = sqrt(2*qe*Te2/me);
% vd1 = veth1*0.0; vd1 = sqrt(2*qe*300/me);
% vd2 = veth1*S; vd2 = sqrt(2*qe*1500/me);
% 
% 
% wpe = sqrt((ne1+ne2)*qe^2/(me*eps)); % Hz
% wpe1 = sqrt(ne1*qe^2/(me*eps)); % Hz
% wpe2 = sqrt(ne2*qe^2/(me*eps)); % Hz
% wpi = sqrt(ni*qe^2/(mi*eps)); % Hz
% Ld = sqrt(veth1^2)/wpe/sqrt(2);
%
% Dispersion relation saved in structure dispcurv


imax = find(dispcurv.wi==max(dispcurv.wi));
vmax = dispcurv.wr(imax)/dispcurv.kvec(imax);


hca = subplot(1,1,1); hold(hca,'on')
lines = plot(hca,dispcurv.kvec*Ld,[dispcurv.wr;dispcurv.wi]);
lines_max = plot(dispcurv.kvec(imax)*Ld,[dispcurv.wr(imax);dispcurv.wi(imax)],'o');
hold(hca,'on')
lines_max(1).Color = lines(1).Color;
lines_max(2).Color = lines(2).Color;
%%
load /Users/Cecilia/MATLAB/cn-matlab/+thesis/dispcurv_modbun.mat
dispcurv1 = dispcurv;
load /Users/Cecilia/MATLAB/cn-matlab/+thesis/dispcurv_ee_s07.mat
dispcurv1 = dispcurv;

wi1 = dispcurv1.wi;
wr1 = dispcurv1.wr;
k1 = dispcurv1.kvec;
vph1 = wr1./k1;
wimax1 = find(wi1==max(wi1));
Ld1 = dispcurv1.Ld;

wi2 = dispcurv2.wi;
wr2 = dispcurv2.wr;
k2 = dispcurv2.kvec;
vph2 = wr2./k2;
wimax2 = find(wi2==max(wi2));
Ld2 = dispcurv2.Ld;



hca = subplot(1,1,1); hold(hca,'on');
lines1 = plot(hca,k1*Ld1,[wr1;wi1]/dispcurv.wpi);
lines_max1 = plot(hca,k1(wimax1)*Ld1,[wr1(wimax);wi1(wimax1)]/dispcurv.wpi,'o');
plot(hca,k1*Ld1,k1*Ld1*0,'-','color',[0.8 0.8 0.8])
hold(hca,'off')
lines_max1(1).Color = lines1(1).Color;
lines_max1(2).Color = lines1(2).Color;

hca.Box = 'on';
hca.FontSize = 14;
hca.YLabel.String = '\omega/\omega_{pi}';
hca.XLabel.String = 'k\lambda_{De}';
text(0.05*max(hca.XLim)/diff(hca.XLim),0.9*max(hca.YLim),{['n_{beam} = ' num2str(dispcurv1.R,'%.2f') ' n_{tot}'],...
            ['v_{beam} = ' num2str(dispcurv1.S,'%.2f') ' v_{te,bg}'],...
            ['v_{te,beam} = ' num2str(dispcurv1.veth2/dispcurv1.veth1,'%.2f') ' v_{te,bg}'],...
            ['v_{ph} = ' num2str(vph1(wimax1)/dispcurv1.veth1,'%.2f') ' v_{te,bg}'],...
            },'fontsize',14,'horizontalalignment','left','verticalalignment','top')
%irf_legend(hca,{'\omega_r'},[0.7,0.6],'color',lines1(1).Color,'fontsize',14)
%irf_legend(hca,{'\gamma'},[0.7,0.3],'color',lines1(2).Color,'fontsize',14)
irf_legend(hca,{'\omega_r'},[0.7,0.75],'color',lines1(1).Color,'fontsize',14)
irf_legend(hca,{'\gamma'},[0.7,0.4],'color',lines1(2).Color,'fontsize',14)


fig = gcf;
fig.Position(4) = fig.Position(4)*0.7;
fig.Position(3) = fig.Position(3)*0.8;

set(gcf, 'InvertHardCopy', 'off');
set(gcf,'paperpositionmode','auto');
set(gcf,'color','white');
%%
hca = subplot(1,1,1); hold(hca,'on');
lines2 = plot(hca,k2*Ld2,[wr2;wi2]/dispcurv.wpi);
lines_max2 = plot(hca,k2(wimax2)*Ld2,[wr2(wimax2);wi2(wimax2)]/dispcurv.wpi,'o');
hold(hca,'on')
lines_max2(1).Color = lines1(1).Color;
lines_max2(2).Color = lines1(2).Color;