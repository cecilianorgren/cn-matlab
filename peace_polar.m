%%
what='DEFlux'
for k=3:4
ic=k;
tin=[cn_toepoch(t1) cn_toepoch(t2)];
eval(['var_name=[''C'',num2str(k),''__CP__PEA__3DXPH__'',what];']);
res=cn_c_caa_construct_subspin_res_data(irf_ssub(['Data__C?_CP_PEA_3DXPH_',what],ic));
%res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXH_DPFlux',ic));
[~,ind]=irf_tlim(res.tt,tin);
specrec=struct('t',res.tt,'dt',res.dtsampling/2,'p_label',['Log ' what ' [' res.dataunits ']']);
specrec.f=res.theta;
specrec.f_label='Pitch angle';
specrec.p=res.pitch_angle(ind,:);
specrec.en=res.en(:);
specrec.data=res.data(ind,:,:);
%%
thi=1:length(res.theta);
r=double(specrec.en);
rlog = log10( r ); % energy levels in eV
r_man = rlog-rlog(1);
r0=rlog(1)-r_man(1);
r_man=[r_man; 2*r_man(end)-r_man(end-1)];
theta = double(specrec.f(thi))+90-(specrec.f(2)-specrec.f(1))/2; % pitch angles
theta=[theta theta(end)+(specrec.f(2)-specrec.f(1))]; % Kolla vinklarna idag
X = r_man*cosd(theta); % xy-coordinates
Y = r_man*sind(theta); % xy-coordinates

resultat=log10(specrec.data(:,thi,:));
%resultat(isnan(resultat))=0;

C=nanmean(resultat,1);
C=nanmean(resultat(:,:,:),1);
%C=resultat(1,:,:);
CC=squeeze(C);

CC=CC';

Xplot=[-flipdim(X(2:end,:),1); X];
Yplot=[flipdim(Y(2:end,:),1); Y];
Cplot=[flipdim(CC(1:end,:),1); CC(1:end,:)];

figure(k+2);
h=surf(Xplot,Yplot,Xplot*0,Cplot(1:end,:));
view(2)
axis equal tight
shading flat 
grid off
cb=colorbar;
ylabel(cb,specrec.p_label)
xlabel('Energy  [keV]')
ylabel('Energy  [keV]')
title(['C',num2str(ic)])

% Ticks
xticks=get(gca,'XTick')-r0;
xticks=[xticks(find(xticks>0)) xticks(end)+1 xticks(end)+2];
xticklabels=cell(size(xticks));
xticklabels={'0.1','1','10'};
xticks=[-flipdim(xticks,2) xticks];
xticklabels=[flipdim(xticklabels,2) xticklabels];
yticks=xticks;
yticklabels=xticklabels;
%xticklabels=num2strget('XTickLabel')-
t1str=datestr(epoch2date(tin(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
t2str=datestr(epoch2date(tin(2)),'MM:SS.FFF');
title([var_name,'    ',t1str,'-',t2str,'UT'])
set(gcf,'PaperPositionMode','auto');
if 0
tick0=log10(ceil(10^r0))-r0;
log10([]);
logscale=log10([2 3 4 5 6 7 8 9 10]);
xtickspos=[];
xticklabels=[];
for k=0:6
    xtickspos=[xtickspos logscale+k*log10(10^k)];
    lstr=eval('10^',num2str(k))
    xticklabels=[ xticklabels; ['', '', '', '', '', '', '', '', '', lstr]'];
    
end

%xtickspos=[logscale, logscale+1*log10(10^k), logscale+2*log10(100), logscale+3*log10(1000), logscale+4*log10(10000)];
xtickspos=xtickspos-tick0;
xticks=[-flipdim(xtickspos,2) xtickspos];
xticklabels=[-flipdim(xticklabels,1) xticklabels];
end
set(gcf,'PaperPositionMode','auto');
set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in','XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)
t1str_p=datestr(epoch2date(tin(1)),'yyyymmddTHHMMSSFFF');
t2str_p=datestr(epoch2date(tin(2)),'MMSSFFF');
eval(['print -depsc2 /Users/Cecilia/Dropbox/Cecilia/EH/',t1str_p,'-',t2str_p,'_PitchDistr',what,'_C',num2str(ic),'.eps']);
end