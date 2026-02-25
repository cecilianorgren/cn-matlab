%% Making figure

% distribution functions
fe = @(v,T,n,vd) cn.maxwellian(v,T(1),n(1),vd(1),'e','1D')+cn.maxwellian(v,T(2),n(2),vd(2),'e','1D'); % electrons
fi = @(v,T,n,vd) cn.maxwellian(v,T,n,vd,'p','1D'); % ions

% distributin parameters
Teb = 60; Tebg = 2000;
R = 0.10;
S = [1.3 0.4];
n = 1;

toPlot = 1;
S = S(toPlot);


% thermal velocities
vtebg = cn_eV2v(Tebg,'eV'); 
vnorm=vtebg;
vlim = 2;
v = linspace(-vlim*vnorm,vlim*vnorm,5000);

% input parameters
Te = [Tebg Teb]; Ti = Tebg;
ne = n*[(1-R) R]; ni = n;
vde = [0 S]*vtebg; vdi = 0;

% phase velocities
vph = [1 0.02]; % in terms of vtebg or S
vph = vph(toPlot); % the case to plot

% trapping velocity
vT = S-vph;
vTvec = linspace(vph-vT,vph+vT,1000);

% plot phase velocities and trapping range
close
set(gcf,'position',[492   585   653   293])

hh=patch([vTvec vTvec([end 1])],[1*fe(vTvec*vnorm,Te,ne,vde)' 0 0]',[1 0.9 0.4],'linewidth',2,'facealpha',1);
hold on;
plot(vph(1)*[1 1],[0 fe(vph(1)*vnorm,Te,ne,vde)],'--k','linewidth',2) 

% plot distributions
finorm = 30;
plot(v/vnorm,fe(v,Te,ne,vde),'k','linewidth',2); 
plot(v/vnorm,fi(v,Ti,ni,vdi)/finorm,'color',[0 0 0],'linewidth',2) 
plot([v([1 end])/vnorm],0*[1 1],'k','linewidth',2) % bottom line


% plot velocity arrows and labels
length = 15; % length of arrow head
fsize = 20; % font size

% label a) and b)
ab_labels = {'a)','b)'};
text(-vlim(1),2*fe(0,Te,ne,vde),ab_labels{toPlot},'fontsize',26) 

% xlabel 'v' with arrow at end of x axis
arrow([vlim*0.8 0.0001],[vlim*1.05 0.0001],'length',15,'baseangle',90,'tipangle',30)
text(vlim*1.07,0,'v_{||}','horizontalalignment','center','verticalalignment','top','fontsize',fsize) 

% label distributions
text(-1,0.8*fe(0,Te,ne,vde),'f_e','fontsize',fsize) 
text(-0.25,1.3*fe(0,Te,ne,vde),'f_i','fontsize',fsize) 




% phase velocity 
if 0 % above
    arrow([vph fe(vph(1)*vnorm,Te,ne,vde)*1.5],[vph fe(vph(1)*vnorm,Te,ne,vde)*1],'length',length)
    text(vph,fe(vph(1)*vnorm,Te,ne,vde)*1.5,'v_{ph}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
else % below    
    arrow([vph fe(0*vnorm,Te,ne,vde)*(-0.2)],[vph fe(0*vnorm,Te,ne,vde)*0],'length',length)
    text(vph,fe(0*vnorm,Te,ne,vde)*(-0.2),'v_{ph}','horizontalalignment','center','verticalalignment','top','fontsize',fsize)
end

% thermal velocity 
% put it at slightly lower than 1, so that it doesn't coincide with vph
vtenorm = 0.9;
if 0 % above
    arrow([vtenorm fe(vph(1)*vnorm,Te,ne,vde)*1.5],[vph fe(vph(1)*vnorm,Te,ne,vde)*1],'length',length)
    text(vtenorm,fe(vph(1)*vnorm,Te,ne,vde)*1.5,'v_{ph}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
elseif 0 % below    
    arrow([vtenorm fe(0*vnorm,Te,ne,vde)*(-0.2)],[vtenorm fe(0*vnorm,Te,ne,vde)*0],'length',length)
    text(vtenorm,fe(0*vnorm,Te,ne,vde)*(-0.2),'v_{te,bg}','horizontalalignment','center','verticalalignment','top','fontsize',fsize)
end

% trapping range
if toPlot == 1; fheight = 1.3; else fheight = 1.8; end
if 0 % above horizaontal arrows    
    arrow([vTvec(1) fheight*fe(0*vnorm,Te,ne,vde)],[vTvec(end) fheight*fe(0*vnorm,Te,ne,vde)],'length',length)
    arrow([vTvec(end) fheight*fe(0*vnorm,Te,ne,vde)],[vTvec(1) fheight*fe(0*vnorm,Te,ne,vde)],'length',length)
    text(vph,fheight*fe(0*vnorm,Te,ne,vde),'2v_{T}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
elseif 0 % below vertical arrows    
    arrow([vTvec(end) fe(vph(1)*vnorm,Te,ne,vde)*(-0.2)],[vTvec(end) fe(vph(1)*vnorm,Te,ne,vde)*0],'length',length)
    text(vTvec(end),fe(vph(1)*vnorm,Te,ne,vde)*(-0.2),'v_{ph}+v_{T}','horizontalalignment','center','verticalalignment','top','fontsize',fsize)
    arrow([vTvec(1) fe(vph(1)*vnorm,Te,ne,vde)*(-0.2)],[vTvec(1) fe(vph(1)*vnorm,Te,ne,vde)*0],'length',length)
    text(vTvec(1),fe(vph(1)*vnorm,Te,ne,vde)*(-0.2),'v_{ph}-v_{T}','horizontalalignment','center','verticalalignment','top','fontsize',fsize)        
else % above vertical arrows    
    arr_start = [vTvec(end) fe(0*vnorm,Te,ne,vde)*fheight];
    arr_end = [vTvec(end) fe(vTvec(end)*vnorm,Te,ne,vde)*1];
    arrow(arr_start,arr_end,'length',length)
    text(arr_start(1),arr_start(2),' v_{ph}+v_{T}','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
    arr_start = [vTvec(1) fe(0*vnorm,Te,ne,vde)*fheight];
    arr_end = [vTvec(1) fe(vTvec(1)*vnorm,Te,ne,vde)*1];
    arrow(arr_start,arr_end,'length',length)
    text(arr_start(1),arr_start(2),'v_{ph}-v_{T}  ','horizontalalignment','center','verticalalignment','bottom','fontsize',fsize)
end    

%plot([vTvec vTvec([end 1])],[fe(vTvec*vnorm,Te,ne,vde)' 0 0]')


axis off
hold off;



