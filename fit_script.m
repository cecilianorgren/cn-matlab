res=c_caa_construct_subspin_res_data(['Data__C4_CP_PEA_3DXPH_PSD']);
t1=[2007 08 31 10 18 40];t2=[2007 08 31 10 18 45];
tint=[toepoch(t1) toepoch(t2)];
[~,ind]=irf_tlim(res.tt,tint);
%%
specrec.f=res.theta;
specrec.f_label='Pitch angle';
specrec.p=res.pitch_angle(ind,:);
specrec.en=res.en(:);
specrec.data=res.data(ind,:,:);
%specrec.bg=bg.data(ind,:,:);
%specrec.bg_en=bg.en(:);
specrec.p_label=['Log PSD [' res.dataunits ']'];                   

specrec.data=squeeze(nanmean(specrec.data,1));
%specrec.bg=squeeze(nanmean(specrec.bg,1));
energy=specrec.en;

%%

if 1 % Select data to plot
    PA0=nanmean(specrec.data(1,:),1)';
    if mod(length(specrec.f),2)
        PA90=nanmean(specrec.data(find(length(specrec.f)),:),1)';
    else
        PA90=nanmean(specrec.data([length(specrec.f)/2 length(specrec.f)/2+1],:),1)';
    end
    PA180=nanmean(specrec.data(end,:),1)';
    %PAbg=nanmean(specrec.bg(:,:),1)'; % One common zero-count level for all levels
end

%
distr=4;
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
    w_n =  [0.006 0.07 0.008 0.001 0.00005]*1e6;        % density m^3
    w_t =  [.07 2 0.35 0.1 0.01];            % temperature [keV]
    w_vd = [0 0 -1 -2.7 -5.3];         % V_drift/V_term
    w_a1 = [0.76 1 1 1 1];           % T_perp/T_par

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
    cmp=[1:2 4 5];
case 3
    % Second fit, slight anisotropy
    w_n =  [0.01 0.05 0.0015 0.003 0.0001]*1e6;        % density m^3
    w_t =  [.04 1.7 0.17 0.1 0.01];            % temperature [keV]
    w_vd = [0 0 -2.6 -1.8 -3.8];         % V_drift/V_term
    w_a1 = [0.76 1 1 1 1];           % T_perp/T_par

    w_m=[0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    cmp=[1 2 3];
case 4
    % Tail, no/bg/tail/drift
    w_n =  [0.01 0.06 0.001 0.00005 0.0001]*1e6;        % density m^3
    w_t =  [.04 2.2 1.4 1.1 0.01];            % temperature [keV]
    w_vd = [0 0 -2.2 -4.8 -3.8];         % V_drift/V_term
    w_a1 = [1 1 1 2 1];           % T_perp/T_par

    w_m=[0 0 0 0 0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
    w_d=[1 1 1 1 1];            % loss cone parameter, 1.0 = no loss cone)
    w_a2=[0 0 0 0 0 0 0];           % ??? (use 0)
    w_pa=[0 90 180];    % pitch angeles
    plotoption=[1];
    titleoption=[0];
    cmp=[ 2 3 ];
        
end

figure(k);
h1=subplot('Position',[0.13 0.288095238095238 0.775 0.636904761904763]);
[h,f,Etot]=whamp.plot_f(w_n(cmp),w_m(cmp),w_t(cmp),w_vd(cmp),...
                            w_d(cmp),w_a1(cmp),w_a2(cmp),w_pa,...
                            plotoption);
                        title(gca,' ')
               
 hold(h1,'off')                    

if 1 % Plotting data
    hold on;
    loglog(h1,specrec.en,PA0,'b*'); %hold('on');
    loglog(h1,specrec.en,PA90,'g*'); %hold(ax,'on');
    loglog(h1,specrec.en,PA180,'r*'); %hold(ax,'on');
    %loglog(h1,specrec.en,PAbg,'k--'); %hold(ax,'on');
    grid(h1,'on')
    set(gca,'xlim',[specrec.en(1)*0.8 specrec.en(end)*1.2])
    set(gca,'ylim',[1e-5 1e3])
    legend('0','90','180')
    %plot(gca,[40 40],[1e-5 1e3])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
    t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
    title(gca,['C',sc,'      ',t1str,' - ',t2str])
    set(gcf,'PaperPositionMode','auto');
    
    hold(h1,'off');
end
 

h2=subplot('Position',[0.13 0.11 0.775 0.122186732186732]);
irf_plot(h2,{cn_toepoch(t1,t2,gsmE3(:,1:2)),cn_toepoch(t1,t2,gsmE4(:,1:2))},'comp')
set(h2,'ColorOrder',[[0 1 0];[0 0 1]]);
irf_legend(gca,{'C3','C4'},[0.02, 0.1]);
ylabel(h2,'E_X GSM')
irf_zoom(h2,'x',tint)