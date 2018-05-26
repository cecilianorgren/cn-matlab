%% Making figure, electron weak beam (bump-on-tail)
Teb = 60;
Tebg = 2000;
vtebg = cn_eV2v(Tebg,'eV');
vteb = cn_eV2v(Teb,'eV');
vnorm = vtebg;

n = 0.06;
v = linspace(-2*vnorm,2*vnorm,5000);
R = 0.10;
S = 1.4;
vph = 1; % in terms of vtebg or S


fbg = cn.maxwellian(v,T,n*(1-R),0,'e','1D'); % background electrons
fb = cn.maxwellian(v,Tb,n*R,vtebg*S,'e','1D'); % fast beam
vnorm=1/vnorm;
%ftot = cn.maxwellian(v,[Teb Tebg],n*[R R-1],[0 vtebg*S],'e','1D'); % fast beam
ftot = fbg+fb;


plot(v*vnorm,fbg+fb,'k','linewidth',2); hold on;
plot([v([1 end])*vnorm],0*[1 1],'k','linewidth',2) % zero line
vbeam = S*vtebg;
vph = 0.5*vbeam;
vtrap = vbeam - vph;

traprange = vph + vtrap*[-1 1];
traprange_vec = [linspace(traprange(1),traprange(2),100) traprange(1)];
ftraprange_vec = [traprange_vec(1:end-1) ];

nottoolow = find(v>traprange(1));
nottoohigh = find(v<traprange(2));
traprange_ind = intersect(nottoolow,nottoohigh);
ftrap = ftot(traprange_ind);

p1 = [v(traprange_ind(1)) v(traprange_ind) v(traprange_ind(end))]'*vnorm;
p2 = [0; ftot(traprange_ind); 0];
p3 = 'b'; repmat([0 0 1],numel(p1),1);
patch(p1,p2,p3)
box off
axis off
%%
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
