%% Making figure, electron weak beam (bump-on-tail)
Teb = 60;
Tebg = 2000;
vtebg = cn_eV2v(Tebg,'eV');
vteb = cn_eV2v(Teb,'eV');
vnorm=vtebg;

n=0.06;
v = linspace(-2*vnorm,2*vnorm,5000);
R = 0.10;
S=  0.4;
vph = 1; % in terms of vtebg or S



fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb = cn.maxwellian(v,Tb,n*R,vtebg*S,'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions
vnorm=1/vnorm;

finorm=1/30;
plot(v*vnorm,fbg+fb,'k','linewidth',2); hold on;

%plot(v*vnorm,fb2,'b','linewidth',2) 
%plot(v*vnorm,fbg,'k','linewidth',2) 
%plot(v*vnorm,fi*finorm,'color',[0 0.7 0],'linewidth',2) 
plot(v*vnorm,fi*finorm,'color',[0 0 0],'linewidth',2) 
plot([v([1 end])*vnorm],0*[1 1],'k','linewidth',2)

fbg0 = cn.maxwellian(0,T,n*(1-R),0,'e','1D');
fbS = cn.maxwellian(S/vnorm,Tb,n*R,vtebg*S,'e','1D'); % fast beam

% phase velocity
fvph = cn.maxwellian(1/vnorm,T,n*(1-R),0,'e','1D') + cn.maxwellian(1/vnorm,Tb,n*R,vtebg*S,'e','1D');
%plot(vph,fvph,'ro') 
plot([vph vph],[fvph fbg0*1.1],'--k') 
plot([S S],[fbS*1.1 fbg0*1.1],'--k') 

% velocity labels

length = 12;
arrow([0 fbg0*1.1],[1 fbg0*1.1],'length',length)
arrow([1 fbg0*1.1],[0 fbg0*1.1],'length',length)
arrow([1.3 fbg0*1.1],[1 fbg0*1.1],'length',length)
arrow([1 fbg0*1.1],[1.3 fbg0*1.1],'length',length)
text(.5,fbg0*1.15,'v_{ph}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
text(1.15,fbg0*1.15,'v_{T}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)

fsize = 20;
%text(S-0.5,max(fb)*1.9,'f_{e,beam}','horizontalalignment','center','fontsize',fsize)
text(-0.7,max(fbg)*0.8,'f_{e}','horizontalalignment','right','fontsize',fsize)
text(-0.1,max(fi)*finorm*0.9,'f_i','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])

axis off
%% Just testing
T = 600;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,1000);
R = 0.15;
S = 1;
vd = vtebg*S;
fbg = cn.maxwellian(v,T,n*(1-R),0,'e');
fb = cn.maxwellian(v,Tb,n*R,vd,'e');
vnorm=1/vtebg;
plot(v*vnorm,fbg,v*vnorm,fb,v*vnorm,fbg+fb)
legend('bg','beam','bg+beam','location','best')
xlabel('v [10^3 km/s]')
set(gca,'xlim',[v(1)*vnorm v(end)*vnorm])

%% Making figure, modified buneman and bump-on-tail, including ions
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,5000);
R = 0.15;
S= [0.4 1.5];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb1 = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D'); % slow beam
fb2 = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions

vnorm=1/vtebg;

plot(v*vnorm,fbg+fb1,'r','linewidth',2); hold on;
plot(v*vnorm,fbg+fb2,'b','linewidth',2) 
plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi/10,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(-0.1,max(fi)/10*0.7,'f_i/10','horizontalalignment','right','fontsize',fsize)
text(S(1)+0.15,max(fb1+fbg)*0.9,'f_{e,beam1}','horizontalalignment','left','fontsize',fsize)
text(S(2)+0.15,max(fb2+fbg)*0.9,'f_{e,beam2}','horizontalalignment','left','fontsize',fsize)
text(-1,max(fbg)*0.8,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
%% Making figure, classical Buneman
Ti = 2000;
Te=Ti/1000;
T = 2000;

vnorm = cn_eV2v(T,'eV');
vtebg = cn_eV2v(Te,'eV');

n=0.06;
v = linspace(-3*vnorm,3*vnorm,5000);
R = 0.15;
%S= [0.4 1.5];

S=1;
fe = cn.maxwellian(v,Te,n,S*vnorm,'e','1D'); % background electrons
fi = cn.maxwellian(v,Ti,n,0,'p','1D'); % ions

vnorm=1/vnorm;
finorm=1/1;
plot(v*vnorm,fe,'b','linewidth',2); hold on;
plot(v*vnorm,fi*finorm,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(-0.1,max(fi)*finorm*0.7,'f_i','horizontalalignment','right','fontsize',fsize)
text(S+0.1,max(fe)*0.8,'f_{e}','horizontalalignment','left','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
%% Making figure, electron two-stream
Ti = 2000;
Te=Ti/1000;
T = 2000;

vnorm = cn_eV2v(T,'eV');
vtebg = cn_eV2v(Te,'eV');

n=0.06;
v = linspace(-3*vnorm,3*vnorm,5000);
R = 0.5;
%S= [0.4 1.5];

S=1;
fe1 = cn.maxwellian(v,Te,n,-S*vnorm,'e','1D'); % background electrons
fe2 = cn.maxwellian(v,Te,n,S*vnorm,'e','1D'); % background electrons
fi = cn.maxwellian(v,Ti,n,0,'p','1D'); % ions

vnorm=1/vnorm;
finorm=1/1;
plot(v*vnorm,fe1,'b','linewidth',2); hold on;
plot(v*vnorm,fe2,'b','linewidth',2);
plot(v*vnorm,fi*finorm,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(-0.1,max(fi)*finorm*0.7,'f_i','horizontalalignment','right','fontsize',fsize)
text(S+0.1,max(fe)*0.8,'f_{e2}','horizontalalignment','left','fontsize',fsize)
text(-S+0.1,max(fe)*0.8,'f_{e1}','horizontalalignment','left','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
%% Making figure, electron weak beam
T = 1200;
Tb = 60;
Tnorm = 2000;
vtebg = cn_eV2v(T,'eV');
vnorm = cn_eV2v(Tnorm,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vnorm,3*vnorm,5000);
R = 0.05;
S=  2.6;


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb = cn.maxwellian(v,Tb,n*R,vtebg*S,'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions
vnorm=1/vnorm;

finorm=1/30;
plot(v*vnorm,fbg+fb,'b','linewidth',2); hold on;
%plot(v*vnorm,fb2,'b','linewidth',2) 
%plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi*finorm,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(S-0.5,max(fb)*1.9,'f_{e,beam}','horizontalalignment','center','fontsize',fsize)
text(-0.7,max(fbg)*0.8,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
text(-0.1,max(fi)*finorm*0.9,'f_i','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
%% Making figure
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,500);
R = 0.15;
S= [0.4 1.5];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D');
fb = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D');
vnorm=1/vtebg;
plot(v*vnorm,fbg+fb,'r','linewidth',2); 
hold on

fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D');
fb = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D');
fi = cn.maxwellian(v,T,n,0,'p','1D');
%plot(v*vnorm,fbg+fb,'color',[0.4 0.4 0.4],'linewidth',2) 
plot(v*vnorm,fbg+fb,'b','linewidth',2) 

plot(v*vnorm,fbg,'k','linewidth',2) 
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
%% Making figure, modified buneman, including ions
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,5000);
R = 0.15;
S= [0.4 1.5];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb1 = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D'); % slow beam
fb2 = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions

vnorm=1/vtebg;

plot(v*vnorm,fbg+fb1,'r','linewidth',2); hold on;
%plot(v*vnorm,fbg+fb2,'b','linewidth',2) 
plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi/10,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(-0.1,max(fi)/10*0.7,'f_i/10','horizontalalignment','right','fontsize',fsize)
text(S(1)+0.15,max(fb1+fbg)*0.9,'f_{e,beam}','horizontalalignment','left','fontsize',fsize)
%text(S(2)+0.15,max(fb2+fbg)*0.9,'f_{e,beam2}','horizontalalignment','left','fontsize',fsize)
text(-1,max(fbg)*0.8,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])

%% Making figure, bump-on-tail, including ions
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,5000);
R = 0.15;
S= [0.4 1.5];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb1 = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D'); % slow beam
fb2 = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions

vnorm=1/vtebg;
%subplot(2,1,1)
%plot(v*vnorm,fbg+fb1,'r','linewidth',2); hold on;
plot(v*vnorm,fbg+fb2,'b','linewidth',2) ; hold on;
plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi/10,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;
text(-0.1,max(fi)/10*0.7,'f_i/10','horizontalalignment','right','fontsize',fsize)
%text(S(1)+0.15,max(fb1+fbg)*0.9,'f_{e,beam1}','horizontalalignment','left','fontsize',fsize)
text(S(2)+0.15,max(fb2+fbg)*0.9,'f_{e,beam}','horizontalalignment','left','fontsize',fsize)
text(-1,max(fbg)*0.8,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])

%% Making figure, mod buneman + bump-on-tail, including ions
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,5000);
R = 0.15;
S= [0.5 1.3];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb1 = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D'); % slow beam
fb2 = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions

vnorm=1/vtebg;
finorm = 20;
%subplot(2,1,1)
plot(v*vnorm,fbg+fb1,'r-','linewidth',2); hold on;
plot(v*vnorm,fbg+fb2,'b-','linewidth',2) ; hold on;
plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi/finorm,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;

text(-0.1,max(fi)/finorm*0.8,['f_i/' num2str(finorm)],'horizontalalignment','right','fontsize',fsize)
text(S(1)+0.15,max(fb1+fbg)*0.9,'f_{e,beam1}','horizontalalignment','left','fontsize',fsize)
text(S(2)+0.15,max(fb2+fbg)*0.9,'f_{e,beam2}','horizontalalignment','left','fontsize',fsize)
text(-0.3,max(fbg)*1.2,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
%xlabel('v/v_{te,bg}')
%xlabel('v/v_{norm}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
set(gca,'xlim',[-1 2.5])
set(gca,'xtick',[])
%% Making figure, mod buneman + bump-on-tail, including ions, no text
T = 2000;
Tb = 60;
vtebg = cn_eV2v(T,'eV');
vteb = cn_eV2v(Tb,'eV');
n=0.06;
v = linspace(-3*vtebg,3*vtebg,5000);
R = 0.15;
S= [0.5 1.3];


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb1 = cn.maxwellian(v,Tb,n*R,vtebg*S(1),'e','1D'); % slow beam
fb2 = cn.maxwellian(v,Tb,n*R,vtebg*S(2),'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions

vnorm=1/vtebg;
finorm = 15;
%subplot(2,1,1)
plot(v*vnorm,fbg+fb1,'r-','linewidth',2); hold on;
plot(v*vnorm,fbg+fb2,'b-','linewidth',2) ; hold on;
plot(v*vnorm,fbg,'k','linewidth',2) 
plot(v*vnorm,fi/finorm,'color',[0 0.7 0],'linewidth',2) 

fsize = 20;

%text(-0.1,max(fi)/finorm*0.8,['f_i/' num2str(finorm)],'horizontalalignment','right','fontsize',fsize)
%text(S(1)+0.15,max(fb1+fbg)*0.9,'f_{e,beam1}','horizontalalignment','left','fontsize',fsize)
%text(S(2)+0.15,max(fb2+fbg)*0.9,'f_{e,beam2}','horizontalalignment','left','fontsize',fsize)
%text(-0.3,max(fbg)*1.2,'f_{e,bg}','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
%xlabel('v/v_{te,bg}')
%xlabel('v/v_{norm}')
%ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])
set(gca,'xlim',[-1 2.5])
set(gca,'xtick',[])
%set()
box(gca,'off')
axis(gca,'off')

%% Making figure, classical Buneman, 
Teb = 60;
Tebg = 2000;
vtebg = cn_eV2v(Tebg,'eV');
vteb = cn_eV2v(Teb,'eV');
vnorm=vtebg;

n=0.06;
v = linspace(-1*vnorm,1*vnorm,5000);
R = 1;
S=  1.7*vteb/vtebg;
vph = 1; % in terms of vtebg or S


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb = cn.maxwellian(v,Tb,n*R,vtebg*S,'e','1D'); % fast beam
fi = cn.maxwellian(v,T,n,0,'p','1D'); % ions
vnorm=1/vnorm;

finorm=1/7.5;
plot(v*vnorm,fbg+fb,'k','linewidth',2); hold on;

%plot(v*vnorm,fb2,'b','linewidth',2) 
%plot(v*vnorm,fbg,'k','linewidth',2) 
%plot(v*vnorm,fi*finorm,'color',[0 0.7 0],'linewidth',2) 
plot(v*vnorm,fi*finorm,'color',[0 0 0],'linewidth',2) 
plot([v([1 end])*vnorm],0*[1 1],'k','linewidth',2)

fbg0 = cn.maxwellian(0,T,n*(1-R),0,'e','1D');
fbS = cn.maxwellian(S/vnorm,Tb,n*R,vtebg*S,'e','1D'); % fast beam

% phase velocity
fvph = cn.maxwellian(1/vnorm,T,n*(1-R),0,'e','1D') + cn.maxwellian(1/vnorm,Tb,n*R,vtebg*S,'e','1D');
%plot(vph,fvph,'ro') 
%plot([vph vph],[fvph fbg0*1.1],'--k') 
%plot([S S],[fbS*1.1 fbg0*1.1],'--k') 

% velocity labels

if 0
length = 12;
arrow([0 fbg0*1.1],[1 fbg0*1.1],'length',length)
arrow([1 fbg0*1.1],[0 fbg0*1.1],'length',length)
arrow([1.3 fbg0*1.1],[1 fbg0*1.1],'length',length)
arrow([1 fbg0*1.1],[1.3 fbg0*1.1],'length',length)
text(.5,fbg0*1.15,'v_{ph}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
text(1.15,fbg0*1.15,'v_{T}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
end
fsize = 20;
%text(S-0.5,max(fb)*1.9,'f_{e,beam}','horizontalalignment','center','fontsize',fsize)
text(0.45,max(fi)*finorm*0.8,'f_{e}','horizontalalignment','left','fontsize',fsize)
text(-0.15,max(fi)*finorm*0.8,'f_i','horizontalalignment','right','fontsize',fsize)
hold off
%legend('bg','beam','bg+beam','location','best')
xlabel('v/v_{te,bg}')
ylabel('f(v_{||})')
set(gca,'xlim',[-2.5 2.5],'yscale','lin','xscale','lin',...
    'yticklabel',[])

axis off

