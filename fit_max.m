%% Make Maxwellian fits to measured distribution functions.
if 0
    product='C4_CP_PEA_3DXPH_PSD';
    res4=c_caa_construct_subspin_res_data(['Data__', product]);            
    bg4=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);
end
%%
ev='C'; t1=[2007 08 31 10 17 50.00]; t2=[2007 08 31 10 17 54.60];
ev='B'; t1=[2007 08 31 10 17 43.50]; t2=[2007 08 31 10 17 47.70]; % many holes
ev='C'; t1=[2007 08 31 10 17 38.00]; t2=[2007 08 31 10 17 42.0];
%%

sc='4';
%res=eval(['res',sc]);
%bg=eval(['bg',sc]);
Vsc=eval(['P',sc]);
Vsc=irf_tlim(Vsc,[toepoch(t1) toepoch(t2)]);
Vsc=mean(Vsc(:,2));
distr=3;
%%
tint=[toepoch(t1) toepoch(t2)];
%%
t1=[2007 08 31 10 17 00.00];
t1=toepoch(t1);
t2=t1+4;
pll=pll+1;
%%
for k=pll
%    tint=[t1 t2];
[~,ind]=irf_tlim(pitch4.t,tint);

if 0
specrec.f=res.theta;
specrec.f_label='Pitch angle';
specrec.p=res.pitch_angle(ind,:);
specrec.en=res.en(:)+Vsc*1.2;
specrec.data=res.data(ind,:,:);
specrec.bg=bg.data(ind,:,:);
specrec.bg_en=bg.en(:);
specrec.p_label=['Log PSD [' res.dataunits ']'];   

specrec.data=squeeze(nanmean(specrec.data,1));
specrec.bg=squeeze(nanmean(specrec.bg,1));
energy=specrec.en;
elseif 1
specrec.f=pitch4.f;
specrec.f_label=pitch4.f_label;
specrec.en=pitch4.en(:)+Vsc*1.2*1;
specrec.data=pitch4.data(ind,:,:);
specrec.bg=pitch4.bg(ind,:,:);
specrec.p_label=pitch4.p_label;                   

specrec.data=squeeze(nanmean(specrec.data,1));
specrec.bg=squeeze(nanmean(specrec.bg,1));
energy=specrec.en; 
elseif 1
specrec.f=hia3.f;
specrec.f_label=hia3.f_label;
specrec.en=hia3.en(:)+Vsc*1.2*0;
specrec.data=hia3.data(ind,:,:);
specrec.bg=hiabg3.data(ind,:,:)*1e20*Units.mp^2/1.602/2;
specrec.p_label=hia3.p_label;    
specrec.bg=1e36*Units.mp^2/(1.602^2)/2./specrec.en;

specrec.data=squeeze(nanmean(specrec.data,1));
%specrec.bg=squeeze(nanmean(specrec.bg,1));
%energy=specrec.en;  
else 
specrec.f=cis4.f;
specrec.f_label=cis4.f_label;
specrec.en=cis4.en(:)+Vsc*1.2*0;
specrec.data=cis4.data(ind,:,:);
%specrec.bg=cis4.bg(ind,:,:);
specrec.p_label=cis4.p_label;                   

specrec.data=squeeze(nanmean(specrec.data,1));
%specrec.bg=squeeze(nanmean(specrec.bg,1));
energy=specrec.en;    
end
%

if 1 % Select data to plot
    PA0=nanmean(specrec.data(1,:),1)';
    if mod(length(specrec.f),2)
        PA90=nanmean(specrec.data(find(length(specrec.f)),:),1)';
    else
        PA90=nanmean(specrec.data([length(specrec.f)/2 length(specrec.f)/2+1],:),1)';
    end
    PA180=nanmean(specrec.data(end,:),1)';
    PAbg=specrec.bg;%nanmean(specrec.bg(:,:),1)'; % One common zero-count level for all levels
end

%%
switch distr
case 1
    % First fit, seems ok: cold, varmer drifting, warm
    w_n=  [0.01 0.009 0.009 0.02]*1e6;       % density m^3
    w_t=  [.1 0.7 0.02 4];            % temperature [keV]
    w_vd= [0 1 0.5 0];         % V_drift/V_term
    w_a1= [1 1 1 1];           % T_perp/T_par

    w_m=  [0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=  [1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2= [0 0 0 0];           % ??? (use 0)
    w_pa= [0 90];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    cmp=[1:2 3 ];
case 2
    % Second fit, slight anisotropy
    % cold background: 0.006 cc, 70 eV
    % hot background: 0.06 cc, 2 keV
    % antiparallel kinda cold beam: 0.008 cc, 300 eV, vd = -2.7 ve
    % antiparallel kinda cold beam: 0.001 cc, 100 eV, vd = -5.3 ve
    % add even colder beam to make up for the little hump:
    %                   0.00005 cc, 10 eV, vd = -5.3 ve
    w_n =  [0.007 0.07 0.008 0.001 0.00005]*1e6;        % density m^3
    w_t =  [.07 2 0.35 0.1 0.01];            % temperature [keV]
    w_vd = [0 0 -1 -2.7 -5.3];         % V_drift/V_term
    w_a1 = [0.7 1 1 1 1];           % T_perp/T_par

    w_m=[0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    
    % calculate current from beams
    bi=4:5;
    v_te = sqrt(2*w_t(bi)*1e3*Units.e/Units.me); % m/s
    v_d = w_vd(bi).*v_te;
    j = Units.e*w_n(bi).*v_d;
    jtot=sum(j*1e9);    % nT
    cmp=[1:2 4:5];
case 3
    % Second fit, slight anisotropy
    w_n =  [0.0028 0.075 0.0005 0.007 0.00017 0.0001]*1e6;        % density m^3
    w_t =  [.06 1.5 0.1 0.1 0.005 0.002];            % temperature [keV]
    w_vd = [-3.5 0 2.1 -1.8 -5.2 -5.3];         % V_drift/V_term
    w_a1 = [1 1 1 1 1 1];           % T_perp/T_par

    w_m=[0 0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    
    % calculate current from beams
    if 0
    bi=[1 2 3 4 5 6];
    v_te = sqrt(2*w_t(bi)*1e3*Units.e/Units.me); % m/s
    v_d = w_vd(bi).*v_te;
    j = Units.e*w_n(bi).*v_d;
    jtot=sum(j*1e9);    % nT
    end
    
    cmp=[1 2];
case 4
    % fourth fit, slight anisotropy
    w_n =  [0.005 0.06 0.0012 0.006 0.005]*1e6;        % density m^3
    w_t =  [.1 1.6 0.15 0.03 0.5];            % temperature [keV]
    w_vd = [-1.8 0 1.7 -1.0 -1.2];         % V_drift/V_term
    w_a1 = [0.8 1 1 1 1];           % T_perp/T_par

    w_m=[0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    
    % calculate current from beams
    bi=[1 2 3 4 ];
    v_te = sqrt(2*w_t(bi)*1e3*Units.e/Units.me); % m/s
    v_d = w_vd(bi).*v_te;
    j = Units.e*w_n(bi).*v_d;
    jtot=sum(j*1e9);    % nT ... no > nA/m2 ?
    
    cmp=[1 2 3 4 ];
case 5 % ions
    % fourth fit, slight anisotropy
    w_n =  [0.01 0.08 0.01 0.08 0.005]*1e6;        % density m^3
    w_t =  [.1 4.6 0.02 0.003 0.5];            % temperature [keV]
    w_vd = [-2.2 0 -2.5 0 -1.2];         % V_drift/V_term
    w_a1 = [0.7 1 1 1 1];           % T_perp/T_par

    w_m=[1 1 1 1 1];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    
    % calculate current from beams
    bi=[1 2 3 4 ];
    v_te = sqrt(2*w_t(bi)*1e3*Units.e/Units.me); % m/s
    v_d = w_vd(bi).*v_te;
    j = Units.e*w_n(bi).*v_d;
    jtot=sum(j*1e9);    % nT ... no > nA/m2 ?
    
    cmp=[1 2 3 ];    
        
end

figure(k);
set(gcf,'position',[900   120   600   800])
h1=subplot('Position',[0.13 0.3 0.775 0.55]);
[h,f,Etot]=irf_whamp_plot_f(w_n(cmp),w_m(cmp),w_t(cmp),w_vd(cmp),...
                            w_d(cmp),w_a1(cmp),w_a2(cmp),w_pa,...
                            plotoption);
                        title(gca,' ')
                        set(h1,'fontsize',16)
                        
               
 hold(h1,'off')                       
%ten=find(Etot(1,:)<max(res.en));
%ten=find(Etot(1,ten)>min(res.en));

if 1 % Plotting data
    hold on;
    loglog(h1,specrec.en,PA0,'b*'); %hold('on');
    loglog(h1,specrec.en,PA90,'g*'); %hold(ax,'on');
    loglog(h1,specrec.en,PA180,'r*'); %hold(ax,'on');
    %loglog(h1,specrec.en,PAbg,'k--'); %hold(ax,'on');
    grid(h1,'on')
    set(gca,'xlim',[1e1 3e4])
    set(gca,'ylim',[1e-3 1e2])
    legend('0','90','180')
    %plot(gca,[40 40],[1e-5 1e3])
    set(gca,'xscale','log','fontsize',16)
    set(gca,'yscale','log','fontsize',16)
    t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
    t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
    title(gca,['C',sc,'      ',t1str,' - ',t2str],'fontsize',16)
    set(gcf,'PaperPositionMode','auto');
    
    hold(h1,'off');
end
if 1
clear infostr;
    nf=0;
    for s=1:length(cmp)
    nf=s;infostr(nf,1)={[num2str(nf),': ',...
        'n_e = ',num2str(sum(w_n(cmp(nf))*1e-6)),' cc, ',...
        'T_e = ', num2str(w_t(cmp(nf))*1e3),' eV, ',...
        'v_d/v_e = ',num2str(sum(w_vd(cmp(nf)))),', ',...
        'T_{perp}/T_{par} = ',num2str(sum(w_a1(cmp(nf)))),...
        ]};    
    end
    
    ann=annotation('textbox',[0.14 0.31 0.64 0.1],...
    'BackgroundColor',[1 1 1],'Color',[0 0 0],...
    'String',infostr,'FontSize',12,...           
    'FitBoxToText','off','LineStyle','none',...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top'); 
end
h2=subplot('Position',[0.13 0.11 0.775 0.122186732186732]);
irf_plot(h2,{cn_toepoch(t1,t2,gsmE3(:,1:2)),cn_toepoch(t1,t2,gsmE4(:,1:2))},'comp')
%irf_zoom(h2,'x',[toepoch(t1) toepoch(t2)])
set(h2,'ColorOrder',[[0 1 0];[0 0 1]]);
irf_legend(gca,{'C3','C4'},[0.02, 0.1]);
ylabel(h2,'E_X GSM')
irf_zoom(h2,'x',tint)
eval(['print -dpng /Users/Cecilia/Dropbox/Cecilia/EH/MaxC',sc,'__',num2str(k)]);

%t1=t1+2;
%t2=t2+2;

%if any(k==600*[10 20 30 40 50 60 70 80 90 100 110]);
%    close all
%end
end
%%
if 0
figure;
ax=gca;

if 1 % Plotting data
    %loglog(ax,specrec.en,PA0,'b*'); hold(ax,'on');
    %loglog(ax,specrec.en,PA90,'g*'); hold(ax,'on');
    loglog(ax,specrec.en,PA180,'r*'); hold(ax,'on');
    loglog(ax,specrec.en,PAbg,'k--'); hold(ax,'on');
    %loglog(ax,Etot(1,ten),f(1,ten),'b'); hold(ax,'on');
    %loglog(ax,Etot(2,ten),f(2,ten),'g'); hold(ax,'on');
    loglog(ax,Etot(3,ten),f(3,ten),'r'); hold(ax,'on');
    
    set(ax,'xlim',[specrec.en(1)*0.8 specrec.en(end)*1.2])
    t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
    t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
    %prod_str{1}=[t1str,'-',t2str,'UT'];
    %prod_str{2}=product;
    %prodstr=[t1str,'-',t2str,' UT    ',product];
    %prodstr(strfind(prodstr,'_'))=' ';
    %irf_legend(ax,{'0^o','90^o','180^o','-- Zero count'},[0.94 0.94]) 
    ylabel(ax,specrec.p_label)
    xlabel(ax,'Energy  [eV]')
    %grid(ax,'off');
    %title(ax,prodstr)
end
end