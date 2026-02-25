% check staff noise level
%% load real magnetic field data
cd /Users/Cecilia/Data/BM/20070831/
sclist = 3;
load mBS
tint=[dBS3(1,1) dBS3(end,1)];
c_eval('B?staff=c_coord_trans(''DSC'',''gsm'',dBS?,''cl_id'',?);',sclist);

c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
c_eval('gsmB?fgm=c_coord_trans(''gse'',''GSM'',gseB?fgm,''cl_id'',?);',sclist);
c_eval('gsmB?fgm=irf_resamp(gsmB?fgm,B?staff);',sclist);
c_eval('gsmB?=c_fgm_staff_combine(gsmB?fgm(:,1:4),B?staff,''plot'');',sclist)

%% make spectra for two different electron hole b-fields
fileToLoad = {'Spis/20140417T233648-2007-08-31-wider-less-dense-3512694.mat',...
              'Spis/20140417T233648-2007-08-31-wider-less-dense-3512694.mat',...
              'Spis/20140409T222215_Tao2011_3051082.mat'};
amplification = [1 10 1];  
nHoles = 1;
changeAmp = 0.5 + 1*rand(1,nHoles);
changeSize = 0.5 + 1*rand(1,nHoles);  
space = ceil(10*rand(1,nHoles));

for pp = 1:numel(fileToLoad)
    % load simulation field data
    varToLoad = {'simBr','simBx','simBy','simBz','modBr','modBx','modBy','modBz','fR','veh','fZ','lr','name','tdB'};
    for ii=1:numel(varToLoad)
        load(['/Users/Cecilia/Research/EH/TestParticleSimulation/' fileToLoad{pp}],varToLoad{ii});
    end
    % take out the data
    % B is in r-th-z system
    r = lr/10; % km
    dR = fR(:,1,1)-r;
    ri = find(abs(dR)==min(abs(dR)));
    % take average of the b-field at this radius
    modBz1D = squeeze(mean(modBz(ri,:,:),2));
    simBz1D = squeeze(mean(simBz(ri,:,:),2));
    simZ = squeeze(fZ(1,1,:));
    % make time series
    simT = simZ/veh;
    tModBz = [simT-simT(1) simBz1D*1e9]; % [t(s) B(nT)]
    % Hz, resample data to this frequency
    if findstr('Tao',name); f = 9000; elseif findstr('2007-08-31',name); f = 450; end     
    T = tModBz(end,1)-tModBz(1,1);
    nt = f*T;
    ttt = linspace(tModBz(1,1),tModBz(end,1),nt)';
    tModBz = irf_resamp(tModBz,ttt);
    %  make spectra
    % replicate time series x times with one signal in between for same time
    % take t1 as some time during 2007-08-31 event
    t1 = toepoch([2007 08 31 10 17 40]); % start time sometime during real event
    dt = diff(tModBz(1:2,1));    
    Bz = [tModBz(:,1)+t1 tModBz(:,2)*0];
    for ii=1:nHoles
        newB = [tModBz(:,1)*changeSize(ii)+Bz(end,1)+dt tModBz(:,2)*changeAmp(ii)];
        Bz = [Bz; newB];        
        for ll=1:space(ii)
            newB = [tModBz(:,1)*changeSize(ii)+Bz(end,1)+dt tModBz(:,2)*0];        
            Bz = [Bz; newB];
        end
    end
    Bz = irf_tappl(Bz,['*' num2str(amplification(pp))]);
    T = Bz(end,1)-Bz(1,1);
    nt = f*T;
    dt = T/nt;
    ttt = tocolumn(Bz(1,1):dt:Bz(end,1));
    %linspace(Bz(1,1),Bz(end,1),nt)';
    Bz = irf_resamp(Bz,ttt);
    smoothBz = [Bz(:,1) smooth(Bz(:,2))];
    % plot spectra
    wds = 2.^(1:14);
    wdstop = size(Bz,1);
    window = wds(find(wds<wdstop,1,'last'));
    smoothSpecrec = irf_powerfft(smoothBz,window,cn.f(Bz),80);
    nosmoothSpecrec = irf_powerfft(Bz,window,cn.f(Bz),80);
    smoothSpecrec.p_label = 'nT^2/Hz';
    nosmoothSpecrec.p_label = 'nT^2/Hz';    
    smoothSpecrec.pav = squeeze(nanmean(smoothSpecrec.p{1},1));
    nosmoothSpecrec.pav = squeeze(nanmean(nosmoothSpecrec.p{1},1));
    
    % collected data
    colBz{pp} = Bz;
    smoBz{pp} = smoothBz;
    smoSpec{pp} = smoothSpecrec;
    nosmoSpec{pp} = nosmoothSpecrec;
    if findstr('Tao',name); 
        names{pp} = ['Tao2011 * ' num2str(amplification(pp))]; 
    elseif findstr('2007-08-31',name); 
        names{pp} = ['2007-08-31 * ' num2str(amplification(pp))]; 
    end         
end
%% plot spectrum of real time series
wd=2^14; % window size
psdh=spectrum.welch('Hann',wd,75);
psdestFGM=psd(psdh,gsmB3fgm(:,2),'Fs',cn.f(gsmB3fgm));
psdestSTA=psd(psdh,B3staff(:,2),'Fs',cn.f(B3staff));

%% plot 3 spectra together
plotSmooth = 0;
if plotSmooth
    plotBz = smoBz;
    spec = smoSpec;
    titleStr = 'smoothed waveforms';
else
    plotBz = colBz;
    spec = nosmoSpec;
    titleStr = 'not smoothed waveforms';
end
nPlot = 4;
for k = 1:nPlot; h(k) = subplot(nPlot,1,k); end
set(h(1),'position',[0.13 0.55 0.7750 0.35])
set(h(2),'position',[0.13 0.39 0.7750 0.09])
set(h(3),'position',[0.13 0.25 0.7750 0.09])
set(h(4),'position',[0.13 0.11 0.7750 0.09])
%h=axes; hold(h,'on');
isub = 0;
if 1 % plot spectra
    isub = isub+1; hca = h(isub); hold(hca,'on');
    nleg = 0;
    if 1 % B_fgm
        plot(hca,psdestFGM.Frequencies, psdestFGM.Data,'k');
        nleg=nleg+1; legends{nleg} = 'FGM';
    end
    if 1 % B_staff        
        plot(hca,psdestSTA.Frequencies, psdestSTA.Data,'b');
        nleg=nleg+1; legends{nleg} = 'STAFF';
    end
    if 1 % B_eh
        pp = 2;
        plot(hca,spec{pp}.f, spec{pp}.pav,'r');
        nleg=nleg+1; legends{nleg} = names{pp};
    end
    if 1 % B_eh
        pp = 1;
        plot(hca,spec{pp}.f, spec{pp}.pav,'g');
        nleg=nleg+1; legends{nleg} = names{pp};
    end
    if 1 % B_eh
        pp = 3;
        plot(hca,spec{pp}.f, spec{pp}.pav,'color',[1 0 1]);
        nleg=nleg+1; legends{nleg} = names{pp};
    end
    title(hca,titleStr)
    ylabel(hca,'nT^2/Hz');   xlabel(hca,'Hz');
    set(hca, 'XScale', 'log','YScale', 'log','ylim',[1e-15 1e5],'xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
    legend(hca,legends)
    hold(hca,'off');
    box(hca,'on')
    grid(hca,'off')
end
if 1 % plot b
    isub = isub+1; hca = h(isub);
    pp = 2;
    plot(hca,plotBz{pp}(:,1)-plotBz{pp}(1,1),plotBz{pp}(:,2),'color','r')%; hold(hca,'on')
    %plot(hca,smoBz{pp}(:,1)-smoBz{pp}(1,1),smoBz{pp}(:,2),'color','k')
    ylabel(hca,'B [nT]'); xlabel(hca,'time [s]')
    set(hca,'ylim',[min(plotBz{pp}(:,2)) max(plotBz{pp}(:,2))])
    grid(hca,'off')
end
if 1 % plot b
    isub = isub+1; hca = h(isub);
    pp = 1;
    plot(hca,plotBz{pp}(:,1)-plotBz{pp}(1,1),plotBz{pp}(:,2),'color','g')
    ylabel(hca,'B [nT]'); xlabel(hca,'time [s]')
    set(hca,'ylim',[min(plotBz{pp}(:,2)) max(plotBz{pp}(:,2))])
    grid(hca,'off')
end
if 1 % plot b
    isub = isub+1; hca = h(isub);
    pp = 3;
    plot(hca,plotBz{pp}(:,1)-plotBz{pp}(1,1),plotBz{pp}(:,2),'color',[1 0 1])
    ylabel(hca,'B [nT]'); xlabel(hca,'time [s]')
    set(hca,'ylim',[min(plotBz{pp}(:,2)) max(plotBz{pp}(:,2))])
    grid(hca,'off')
end
    



