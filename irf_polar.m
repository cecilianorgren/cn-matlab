function irf_polar(varargin)

% Plots a polar distribution of particle data.
% First argument is time interval in epoch.
% If two products are given, plot together in 
% opposite side of plot.

time=varargin{1};
switch size(varargin,2)
    case 2
        product=varargin{2};
        if any(strfind(product,'3DXPH'))
            res=cn_c_caa_construct_subspin_res_data(product);
            res
            [~,ind]=irf_tlim(res.tt,time);
            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            specrec.en=res.en(:);
            specrec.data=res.data(ind,:,:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];
            
        elseif any(strfind(product,'3DRL'))


        elseif any(strfind(product,'PITCH_3DXH')) || any(strfind(product,'PITCH_3DRH')) || any(strfind(product,'PITCH_3DRL'))
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',product]);
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
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            
            [~,ind]=irf_tlim(tt,time);
            
            specrec.f_label='Pitch angle';
            specrec.data=specrec.data(ind,:,:);
            
            thi=1:length(theta);
            r=double(specrec.en)';
            rlog = log10( r ); % energy levels in eV
            r_man = (rlog-rlog(end)+0.01);
            r0=rlog(1)-r_man(1);
            theta = double(specrec.f(thi))+90; % pitch angles
            X = r_man*cosd(theta); % xy-coordinates
            Y = r_man*sind(theta); % xy-coordinates

            resultat=log10(specrec.data(:,thi,:));

            C=nanmean(resultat,1);
            CC=double(squeeze(C));
            CC=CC';

            Xplot=[-flipdim(X,1) X];
            Yplot=[flipdim(Y,1) Y];
            Cplot=[flipdim(CC,1) CC];

            figure;
            h=pcolor(Xplot,Yplot,Cplot);
            axis equal tight 
            cb=colorbar;
            %ylabel(cb,specrec.p_label)
            shading flat
            xlabel('Energy [keV]')
            ylabel('Energy [keV]')
            title(product)
            
        elseif any(strfind(product,'3DXPH'))


        end
    case 3
        
end
if 0
%%
what='DEFlux';
prod='3DXPH';
%prod='3DRL';
for k=3:4
ic=k;
tin=[cn_toepoch(t1) cn_toepoch(t2)];
eval(['var_name=[''C'',num2str(k),''__CP__PEA__'',prod,''__'',what];']);
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
r=[r(1)-(r(2)-r(1))/2; r(1:end-1)+(r(2:end)-r(1:end-1))/2 ;r(end)+(r(end)-r(end-1))/2];
rlog = log10( r ); % energy levels in eV
r_man = rlog-rlog(1);
r0=rlog(1)-r_man(1);
%r_man = [r_man(1)]
[r_man; 2*r_man(end)-r_man(end-1)];
theta = double(specrec.f(thi))-90-(specrec.f(2)-specrec.f(1))/2; % pitch angles
theta=[theta theta(end)+(specrec.f(2)-specrec.f(1))]; % Kolla vinklarna idag
X = r_man*cosd(theta); % xy-coordinates
Y = r_man*sind(theta); % xy-coordinates

resultat=log10(specrec.data(:,thi,:));
%resultat(isnan(resultat))=0;

C=nanmean(resultat,1);
%C=nansum(resultat,1);
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
if strcmp(what,'PSD'); caxis(gca,[-4 1]); end

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
end
end

function [tt,dtsampling]=subspintime(dataobject,phi)
% construct subspin time vector
% phi are azimuthal angles (spin period is divided in the number of azimuth
% angles)
timevar=getv(dataobject,dataobject.VariableAttributes.DEPEND_0{1,2});
tt=timevar.data(:);
tt=repmat(tt,1,length(phi));

if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
    if ischar(timevar.DELTA_PLUS)
        deltaplus= getv(dataobject,timevar.DELTA_PLUS);
        dtplus=deltaplus.data(1,:);
        if dtplus>5, % temporary solution for CIS problems
            if (dtplus/2>3.5) && (dtplus/2 < 4.5), dtplus=dtplus/2;
            elseif (dtplus/3>3.5) && (dtplus/3 < 4.5), dtplus=dtplus/3;
            elseif (dtplus/4>3.5) && (dtplus/4 < 4.5), dtplus=dtplus/4;
            end
        end
    elseif isnumeric(timevar.DELTA_PLUS)
        dtplus=timevar.DELTA_PLUS;
    end
    if ischar(timevar.DELTA_MINUS)
        deltaminus= getv(dataobject,timevar.DELTA_MINUS);
        dtminus=deltaplus.data(1,:);
    elseif isnumeric(timevar.DELTA_MINUS)
        dtminus=timevar.DELTA_MINUS;
    end
else
    dtplus=2;
    dtminus=2;
end
spin_period=dtplus+dtminus;
dtsampling=spin_period/length(phi);
for j=length(phi):-1:1,
    tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
end
tt=reshape(tt',numel(tt),1);
end