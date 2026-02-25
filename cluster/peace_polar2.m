%%
for k=3:4
ic=k;
tin=[cn_toepoch(t1) cn_toepoch(t2)];
eval(['var_name=[''C'',num2str(k),''__CP__PEA__PITCH__3DXH__DEFlux''];']);
c_eval('[caaData,dataobject,Data,Data_units]=c_caa_var_get(''Data__C?_CP_PEA_PITCH_3DXH_DEFlux'');',ic)
theta=Data.dep_x{2}.data(1,:);
t=Data.t;
en=Data.dep_x{3}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
dataraw=Data.data;
dataraw(:,:,:,nan_en)=[];
dataraw(:,nan_phi,:,:)=[];

%dataraw(dataraw==caaData.FILLVAL)=NaN;
dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
tt=subspintime(dataobject,phi);
repmat(t(:),1,length(phi));
specrec.f=theta;
specrec.data=data;
specrec.en=en;
specrec.t=tt;
%%
%Data.data
%res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXH_DEFlux',ic));
%res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXH_DPFlux',ic));
[~,ind]=irf_tlim(tt,tin);
%specrec=struct('t',res.tt,'dt',res.dtsampling/2,'p_label',['Log DEFlux [' res.dataunits ']']);
%specrec.f=res.theta;
specrec.f_label='Pitch angle';
%specrec.p=res.pitch_angle(ind,:);
%specrec.en=res.en(:);
specrec.data=specrec.data(ind,:,:);
%%
thi=1:length(theta);
r=double(specrec.en)';
rlog = log10( r ); % energy levels in eV
r_man = (rlog-rlog(end)+0.01);
r0=rlog(1)-r_man(1);
theta = double(specrec.f(thi))-90; % pitch angles
X = r_man*cosd(theta); % xy-coordinates
Y = r_man*sind(theta); % xy-coordinates

resultat=log10(specrec.data(:,thi,:));
%resultat(isnan(resultat))=0;

C=nanmean(resultat,1);
C=nanmean(resultat(:,:,:),1);
%C=resultat(1,:,:);
CC=double(squeeze(C));

CC=CC';

Xplot=[-flipdim(X,1) X];
Yplot=[flipdim(Y,1) Y];
Cplot=[flipdim(CC,1) CC];

figure(k);
h=pcolor(Xplot,Yplot,Cplot);
axis equal tight 
cb=colorbar;
%ylabel(cb,specrec.p_label)
shading flat
xlabel('Energy  [eV]')
ylabel('Energy [keV]')
title(['C',num2str(ic)])

% Ticks
if 0
xticks=get(gca,'XTick')-r0;
xticks=[xticks(find(xticks>0)) xticks(end)+1 xticks(end)+2];
xticklabels=cell(size(xticks));
xticklabels={'0.1','1','10'};
xticks=[-flipdim(xticks,2) xticks];
xticklabels=[flipdim(xticklabels,2) xticklabels];
yticks=xticks;
yticklabels=xticklabels;
end
%xticklabels=num2strget('XTickLabel')-
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
t1str=datestr(epoch2date(tin(1)),'dd-mmm-yyyy__HH:MM:SS.FFF');
t2str=datestr(epoch2date(tin(2)),'MM:SS.FFF');
title([var_name,'    ',t1str,'-',t2str,'UT'])
set(gcf,'PaperPositionMode','auto');
%set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in','XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)
eval(['print -depsc2 /Users/Cecilia/Dropbox/Cecilia/EH/',t1str,'-',t2str,'PitchDistr_C',num2str(ic),'.eps']);
end
