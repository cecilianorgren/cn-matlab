
%% first things first

tint = irf.tint('2017-01-04T05:03:18/2017-01-04T05:03:43');
ic = 2;

% eis mode

eisMode = 'brst';


%% get fgm data

B = mms.get_data('B_gse_fgm_brst_l2',tint,ic);

%% get fpi distribution

iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);

iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
iPDist.data = iPDist.data; % to SI


%% get EIS diff particle flux of protons

% unit of Eeis is [eV] and unit of EISdpfX is []

switch eisMode
    case 'srvy'
        c_eval('EISdpf! = mms.db_get_ts(''mms?_epd-eis_srvy_l2_phxtof'',''mms?_epd_eis_phxtof_proton_P4_flux_t!'',tint);',1:4,0:5)
        
        
        % assume all energy tables are the same        
        [filepath,filename] = mms.get_filepath('mms2_epd-eis_srvy_l2_phxtof',tint(1));
        do = dataobj([filepath,filename]);
        Eeis = 1e3*do.data.mms2_epd_eis_phxtof_proton_t0_energy.data;
        
        %Eeis = 1e3*[12.0502863036527,17.3272271295785,23.7937050498656,31.2732089814909,41.4035650638927,54.3957608642438,70.5562470657276];
        
    case 'brst'
        c_eval('EISdpf! = mms.db_get_ts(''mms?_epd-eis_brst_l2_phxtof'',''mms?_epd_eis_brst_phxtof_proton_P4_flux_t!'',tint);',1:4,0:5)
        %         Eeis = 1e3*[10.5351791631657,11.5243741695591,12.6127950907006,...
        %             13.8078238690968,15.0820555628203,16.4947264299836,18.0496315751873,...
        %             19.7377658310154,21.6108491272324,23.5828327390922,25.7281732736702,...
        %             28.1473266936667,30.7321623059621,33.6260761909563,36.8466769193641,...
        %             40.1890308417608,44.1347074307738,49.8928050925181,56.5827447195698,...
        %             63.4273009042636,71.8381681280919];
        
        [filepath,filename] = mms.get_filepath('mms2_epd-eis_brst_l2_phxtof',tint(1)+20);
        do = dataobj([filepath,filename]);
        Eeis = 1e3*do.data.mms2_epd_eis_brst_phxtof_proton_t0_energy.data;
end



%% Convert EIS diff particle flux to phase-space density


mm = 1; % ion mass
% energy table for each time step and energy index
ETabEis = repmat(Eeis,EISdpf0.length,1);
% copy objects
c_eval('EISpsd? = EISdpf?;',0:5)
% sort of magic conversion :(
c_eval('EISpsd?.data = EISdpf?.data/1e12*mm^2*0.53707./ETabEis;',0:5)
% mysterious correction factor
c_eval('EISpsd?.data = EISpsd?.data*1e-3;',0:5)



%% average EIS psd from all detectors (same time stamps)

EISpsd = (EISpsd0+EISpsd1+EISpsd2+EISpsd3+EISpsd4+EISpsd5)/6;


%% EIS average spectrogram structure

% time step of eis detectors 
dTEis = median(diff(EISdpf0.time.epochUnix));

EISspecrec = [];
EISspecrec.p = EISpsd.data;
EISspecrec.t = EISpsd.time.epochUnix+dTEis/2;
EISspecrec.f = Eeis*1e-3;
EISspecrec.f_label = 'E (keV)';


%% plot overview time series

h = irf_plot(3,'newfigure');


% B 1sc
hca = irf_panel(h,'bxyz');
irf_plot(hca,B)
hold(hca,'on')
irf_plot(hca,B.abs)
ylabel(hca,'$\mathbf{B}$ [nT]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$B_x$';'$B_y$';'$B_z$';'$|B|$'},[1.02,0.98],'Fontsize',15,'interpreter','latex')
irf_legend(hca,['MMS',num2str(ic)],[0.02,1.02],'Fontsize',15,'interpreter','latex');



% fi eis
hca = irf_panel(h,'fieis');
irf_spectrogram(hca,EISspecrec,'donotshowcolorbar')
hca.YScale = 'log';
ylabel(hca,'E [keV]','fontsize',15,'interpreter','latex')
irf_legend(hca,'epd{\textunderscore}eis{\textunderscore}phxtof{\textunderscore}proton{\textunderscore}P4',[0.02,0.98],'Fontsize',15,'interpreter','latex');

% fi fpi
hca = irf_panel(h,'fifpi');
iPDistSI = iPDist;
iPDistSI.data = iPDist.data*1e12;
irf_spectrogram(hca,iPDistSI.omni.specrec,'donotshowcolorbar')
hca.YScale = 'log';
irf_legend(hca,'FPI-DIS',[0.02,0.98],'Fontsize',15,'interpreter','latex');
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];

% match color limits
h(2).CLim = h(3).CLim; pause(0.01)

hcb = colorbar(h(2));
hcb.Position(1) = .85;
ylabel(hcb,'$\log{f_i}$ [s$^3$\,m$^{-6}$]','fontsize',14,'interpreter','latex')
hcb.LineWidth = 1.3;

hcb = colorbar(h(3));
hcb.Position(1) = .85;
ylabel(hcb,'$\log{f_i}$ [s$^3$\,m$^{-6}$]','fontsize',14,'interpreter','latex')
hcb.LineWidth = 1.3;

% more stuff
irf_colormap(h(1),'waterfall')
irf_colormap(h(2),'waterfall')


irf_zoom(h,'x',tint+[-dTEis/2,dTEis/2])
irf_plot_axis_align(h)

pause(0.01)
for jj = 1:length(h)
    irf_zoom(h(jj),'y',h(jj).YLim)
    h(jj).Position(1) = 0.13;
    h(jj).Position(3) = 0.70;
    h(jj).Layer = 'top';
    h(jj).FontSize = 15;
    h(jj).LineWidth = 1.3;
    h(jj).YLabel.Position(1) = -0.12;
end

%% Resample FPI data to EIS time resolution 
% NOTE:
% only for non-alternating energy table of FPI!!!

emat = iPDist.energy; % in eV

FPIspecrec = EISspecrec;
FPIspecrec.f = emat(1,:); % assume constant!

FPIspecrec.p = zeros(EISdpf0.length,size(emat,2)); 

for ii = 1:EISdpf0.length-1
    
    tintTemp = EISdpf0.time(ii:ii+1);
    iPDistSItemp = iPDistSI.tlim(tintTemp);
    fpiOmniTemp = iPDistSItemp.omni;
    
    FPIspecrec.p(ii,:) = nanmean(fpiOmniTemp.data);
    
    
end

% special care for last point

tintTemp = EISdpf0.time(end)+[0,dTEis];
iPDistSItemp = iPDistSI.tlim(tintTemp);
fpiOmniTemp = iPDistSItemp.omni;

FPIspecrec.p(end,:) = nanmean(fpiOmniTemp.data);


%% check FPI averages

h = irf_plot(2,'newfigure');

% fi fpi full
hca = irf_panel(h,'fifpifull');
iPDistSI = iPDist;
iPDistSI.data = iPDist.data*1e12;
irf_spectrogram(hca,iPDistSI.omni.specrec,'donotshowcolorbar')
hca.YScale = 'log';
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];

% fi fpi avg
hca = irf_panel(h,'fifpiavg');
iPDistSI = iPDist;
iPDistSI.data = iPDist.data*1e12;
irf_spectrogram(hca,FPIspecrec,'donotshowcolorbar')
hca.YScale = 'log';
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];

hcb = colorbar(h(1));
hcb.Position(1) = .85;
ylabel(hcb,'$\log{f_i}$ [s$^3$\,m$^{-6}$]','fontsize',14,'interpreter','latex')
hcb.LineWidth = 1.3;

hcb = colorbar(h(2));
hcb.Position(1) = .85;
ylabel(hcb,'$\log{f_i}$ [s$^3$\,m$^{-6}$]','fontsize',14,'interpreter','latex')
hcb.LineWidth = 1.3;

% more stuff
irf_colormap(h(1),'waterfall')
irf_colormap(h(2),'waterfall')

irf_zoom(h,'x',tint+[-dTEis/2,dTEis/2])
irf_plot_axis_align(h)
hcb.LineWidth = 1.3;
pause(0.01)
for jj = 1:length(h)
    irf_zoom(h(jj),'y',h(jj).YLim)
    h(jj).Position(1) = 0.13;
    h(jj).Position(3) = 0.70;
    h(jj).Layer = 'top';
    h(jj).FontSize = 15;
    h(jj).LineWidth = 1.3;
    h(jj).YLabel.Position(1) = -0.12;
end



%% make line plots comparing FPI to EIS for all time steps

cmap = jet;
col = interp1(linspace(1,EISdpf0.length,64),cmap,1:EISdpf0.length);


fig = figure;

hca = axes(fig);
hold(hca,'on')

% for legend
hlTemp1 = plot(hca,FPIspecrec.f,FPIspecrec.p(1,:),'-v','color','k','linewidth',2);
hlTemp2 = plot(hca,EISspecrec.f,EISspecrec.p(1,:),'-o','color','k','linewidth',2);

for ii = 1:EISdpf0.length-1
    plot(hca,FPIspecrec.f,FPIspecrec.p(ii,:),'-v','color',col(ii,:),'linewidth',2)
    plot(hca,EISspecrec.f*1e3,EISspecrec.p(ii,:),'--o','color',col(ii,:),'linewidth',2)
end

hca.XLim(1) = 1e2; 
hca.XScale = 'log';
hca.YScale = 'log';


hleg = legend(hca,'FPI','EIS');
hleg.FontSize = 15;
hleg.Interpreter = 'latex';

%delete(hlTemp1)
%delete(hlTemp2)

hcb = colorbar(hca);
colormap(hca,col);
hca.CLim = [0,diff(EISdpf0.time([1,end]).epochUnix+dTEis)];

ylabel(hca,'$f_i$ [s$^3$\,m$^{-6}$]','fontsize',15,'interpreter','latex')
xlabel(hca,'$E$ [eV]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$T$ [s]','fontsize',15,'interpreter','latex')

hca.Box = 'on';
hca.LineWidth = 1.3;
hcb.LineWidth = 1.3;



%% make line plots comparing FPI to EIS for one time step 

it = 9;

fig = figure;

hca = axes(fig);

%%
hold(hca,'on')

% for legend
%hlTemp1 = plot(hca,FPIspecrec.f,FPIspecrec.p(it,:),'-','color','k','linewidth',2);
%hlTemp2 = plot(hca,EISspecrec.f,EISspecrec.p(it,:),'-o','color','k','linewidth',2);

plot(hca,FPIspecrec.f,FPIspecrec.p(it,:),'-v','color','k','linewidth',2)
plot(hca,EISspecrec.f*1e3,EISspecrec.p(it,:),'--o','color','k','linewidth',2)


hca.XLim(1) = 2e2;
hca.XScale = 'log';
hca.YScale = 'log';


hleg = legend(hca,'FPI','EIS');
hleg.FontSize = 15;
hleg.Interpreter = 'latex';


ylabel(hca,'$f_i$ [s$^3$\,m$^{-6}$]','fontsize',15,'interpreter','latex')
xlabel(hca,'$E$ [eV]','fontsize',15,'interpreter','latex')

hca.Box = 'on';
hca.LineWidth = 1.3;

%%
title(hca,[EISpsd.time(it).toUtc,'+',num2str(dTEis),' s'])

