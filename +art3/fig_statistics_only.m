load('/Users/Cecilia/Research/EH2/wavespeeds.mat')
yES = wavespeeds.vtrESwaves./wavespeeds.vthESwaves;
xES = wavespeeds.vESwaves./wavespeeds.vthESwaves;
yEH = wavespeeds.vtreholes./wavespeeds.vtheholes;
xEH = wavespeeds.veholes./wavespeeds.vtheholes;



% additional good options
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
hold(gca,'on')

hleg(1)=plot(xES,yES,'rx','linewidth',1.5,'markersize',5);
hleg(2)=plot(xEH,yEH,'ko','linewidth',1.5,'markersize',5); 
set(gca,'xlim',[2e-3 2e0],'ylim',[1e-2 3e1],'xtick',[1e-2 1e-1 1e0])
xlabel(gca,'v_{ph}/v_{te,bg}')
ylabel(gca,'v_{T}/v_{te,bg}')
%axis square
%axis equal
set(gca,'yscale','log','xscale','log') 
set(gca,'xlim',[2e-3 1e0],'ylim',[1e-2 1e0],'xtick',[1e-2 1e-1 1e0])
%axis equal
box on
hleg(3) = loglog([1e-3 1e1],[1e-3 1e1],'color',[1 0.7 0.1],'linewidth',1.5);
hlegg = legend([hleg(2) hleg(1) hleg(3)],'Electrostatic solitary waves','Electrostatic waves','v_T=v_{ph}','location','northwest');
hlegg.Box = 'off';
hold(gca,'off')
%%
hleg(4)=loglog(xP,yP,'color',[0.5 0.8 0.1],'linewidth',2); hold(gca,'on')
%hleg(5)=loglog(xPow,yPow,'color',[0.3 0.7 0.9],'linewidth',2); hold(gca,'on')
hleg(6)=loglog(xTheory(:),yTheory(:),'.','color',[0.7 0.7 0.7]);
hleg(7)=loglog(xPthe,yPthe,'color',[1 0.7 0.1],'linewidth',2); hold(gca,'on')
hleg(8)=loglog(xMT,yMT,'^','color',[0.5 0.2 1],'linewidth',2); hold(gca,'on')
hleg(9)=loglog(xEH(40),yEH(40),'^','color',[0.5 0.2 1],'linewidth',2,'markersize',5); 
hold(gca,'off')
%hleg(2)=loglog(xES,yES,'rx','markersize',7);
%hleg(3)=loglog(xEH,yEH,'ko','markersize',7); hold(gca,'off')

set(gca,'xlim',[2e-3 2e0],'ylim',[1e-2 3e1],'xtick',[1e-2 1e-1 1e0])
xlabel(gca,'v/v_{te,bg}')
ylabel(gca,'v/v_{T}')
%legend(hca,strMod,'location','northwest')
hlegg=legend([hleg(2) hleg(3)],'ES','ESW','location','northwest');
set(hlegg,'box','off')

if 0
    hlegg2 = legend(hleg(1),strMod,'location','southeast');
    set(hlegg2,'box','off')
end
%
axis square
set(gca,'yscale','log','xscale','log') 
box on
hold off
%% make vT
% run combine_dg_cn.m first to get base data
vdmat = repmat(s{4}.S,numel(s{4}.R),1)*s{4}.vtebg;
vT = vdmat-s{4}.vph;

% pcolor(s{4}.R,s{4}.S,vT')

plR = find(R>0.1 & R<0.6); 

xTheory = s{4}.vph(plR,:)/s{4}.vtebg;
yTheory = s{4}.vph(plR,:)./vT(plR,:);


% Make polynomial fit of theoretical log values
P = polyfit(log(xTheory),log(yTheory),1);
xP = ppxx;
yP = exp(P(1)*log(xP)+P(2));
fittP = fit(xP',yP','power1');
strP = ['fit: y = ' num2str(fittP.a) '*x^{' num2str(fittP.b) '}'];

hold on;
%loglog(xTheory,yTheory,'.','color',[0.7 0.7 0.7])

%% different comparison, more info
% normalize everyting to vtebg instead
plR = find(s{4}.R>0.0 & s{4}.R<0.53); 
plS =  find(s{4}.S>0.0 & s{4}.S<2); 
yES = wavespeeds.vtrESwaves./wavespeeds.vthESwaves;
xES = wavespeeds.vESwaves./wavespeeds.vthESwaves;
yEH = wavespeeds.vtreholes./wavespeeds.vtheholes;
xEH = wavespeeds.veholes./wavespeeds.vtheholes;
yTheory = vT(plR,plS)/s{4}.vtebg;
xTheory = s{4}.vph(plR,plS)/s{4}.vtebg;

strP = ['x=y5'];


hleg(2)=loglog(xES,yES ,'rx','linewidth',1.5,'markersize',5); hold on
hleg(3)=loglog(xEH,yEH,'ko','linewidth',1.5,'markersize',5); 
hleg(4)=loglog([1e-3 1e1],[1e-3 1e1],'g');
hleg(5)=loglog([1e-3 1e1]*0.2,[1e-3 1e1],'b');
hleg(6)=loglog([1e-4 1e1],1*veth2/veth1*[1 1],'m');
hhh=loglog(xTheory,yTheory,'-','color',[1 0.6 0.6]); hleg(7)=hhh(1);
hhh=loglog(xTheory',yTheory','-','color',[0.6 0.6 1]); hleg(8)=hhh(1);
loglog(xTheory,yTheory,'o','color',[0.7 0.7 0.7]);
hlegg=legend([hleg(2) hleg(3) hleg(4) hleg(5) hleg(6) hleg(7) hleg(8)],'ES','ESW','v_T=v_{ph}','v_T=5v_{ph}','v_{te,beam}/v_{te,bg}','S=const.','R=const.','location','northwest');
set(hlegg,'box','off')

hold off;
ylabel(gca,'v_T/v_{te,bg}')
xlabel(gca,'v_{ph}/v_{te,bg}')
set(gca,'ylim',[1e-2 2e0],'xlim',[3e-3 2e0])
%axis equal
axis square

%% different comparison, less info, for paper
% normalize everyting to vtebg instead
plR = find(s{4}.R>0 & s{4}.R<0.5); 
plS =  find(s{4}.S>0.0 & s{4}.S<2.0); 
yES = wavespeeds.vtrESwaves./wavespeeds.vthESwaves;
xES = wavespeeds.vESwaves./wavespeeds.vthESwaves;
yEH = wavespeeds.vtreholes./wavespeeds.vtheholes;
xEH = wavespeeds.veholes./wavespeeds.vtheholes;
yTheory = vT(plR,plS)/s{4}.vtebg;
xTheory = s{4}.vph(plR,plS)/s{4}.vtebg;

% Magnetotail dots
qe = 1.6022e-19;
me = 9.1090e-31;
phi_MT = [1800 180];
Te_MT = 4000;
v_te_MT = sqrt(2*qe*Te_MT/me);
vT_MT = sqrt(2*qe*phi_MT/me);
v_ph_MT = [14000 1400]*1e3;
xMT = v_ph_MT/v_te_MT;
yMT = v_ph_MT./vT_MT;

strP = ['x=y5'];

hleg(1)=loglog(xEH(40),yEH(40),'^','color',[0.4 0.4 1],'linewidth',2,'markersize',8,'markerfacecolor',[0.4 0.4 1]);  hold on
hleg(2)=loglog(xES,yES ,'rx','linewidth',2,'markersize',7); hold on
hleg(3)=loglog(xEH,yEH,'ko','linewidth',2,'markersize',7); 
hleg(4)=loglog([1e-3 1e1],[1e-3 1e1],'color',[1 0.7 0.1],'linewidth',1.5);
%hleg(5)=loglog([1e-3 1e1]*0.2,[1e-3 1e1],'b','linewidth',1.5);
%hleg(6)=loglog([1e-4 1e1],1*veth2/veth1*[1 1],'m');
%hhh=loglog(xTheory,yTheory,'-','color',[1 0.6 0.6]); hleg(7)=hhh(1);
%hhh=loglog(xTheory',yTheory','-','color',[0.6 0.6 1]); hleg(8)=hhh(1);
hhh=loglog(xTheory,yTheory,'o','color',[0.7 0.7 0.7]);  hleg(8)=hhh(1);
hhh=loglog(xMT,yMT,'^','color',[0.2 0.8 0.2],'linewidth',1.5,'markersize',8,'markerfacecolor',[0.2 0.8 0.2]);  hleg(9)=hhh(1);

hlegg=legend([hleg(8) hleg(2) hleg(3) hleg(1) hleg(9) hleg(4)],'Model','ES','ESW','MP event','MT events','v_T=v_{ph}','location','northeastoutside');
set(hlegg,'box','off')

hold off;
ylabel(gca,'v_T/v_{te,bg}','fontsize',20)
xlabel(gca,'v_{ph}/v_{te,bg}','fontsize',20)
set(gca,'ylim',[1e-2 1e0],'xlim',[3e-3 2e0],'fontsize',16)
%axis equal
axis square





