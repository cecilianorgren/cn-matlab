cd /Users/Cecilia/Data/BM/20070831
t1=[2007 08 31 10 13 00];
t2=[2007 08 31 10 25 00];
tint=toepoch([t1;t2])';

%% Anisotropy data   
fig1 = 1;
if fig1
    loadfrommat=0;
    doDEFlux = 0;
    doDPFlux = 0;
    doPSD = 1;
    
    switch loadfrommat
        case 1
            load an3
        case 0
            sc = 4;
            if doPSD; varname=irf_ssub('C?_CP_PEA_PITCH_3DXH_PSD',sc);
            elseif doDEFlux; varname=irf_ssub('C?_CP_PEA_3DXPH_DEFlux',sc);
            else varname=irf_ssub('C?_CP_PEA_3DXPH_DPFlux',sc); end

            xtra_vn=''; % could be _HEEA or _LEEA when using pitch full
            [or_var,dobj,varmat,varunits]=c_caa_var_get(['Data',xtra_vn,'__',varname]);
            [or_bg,dobj_bg,varmat_bg,varunits_bg]=c_caa_var_get(['BackgroundLevel',xtra_vn,'__',varname]);
            varname2=irf_ssub(['Sweep_Energy',xtra_vn,'__',varname],sc);
            %varname3=irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
            energies=c_caa_var_get(varname2,'mat');
            energy_units=c_caa_var_get(varname2,'units'); % These give eV, should be keV
            %piichangle=c_caa_var_get(varname3,'mat');

            % Taking away all fluxes that are below one, these are something else,
            % what?
            %
            var=or_var;
            var.data(or_var.data == or_var.FILLVAL)=NaN;

            % taking mean of azimuthal angles
            deflux=squeeze(nanmean(var.data(:,:,:,:),2));

            % Making the basic flux structure
            b_flux.t=energies(:,1);
            b_flux.f=energies(1,2:end); nan_en=~isnan(b_flux.f);
            b_flux.f=b_flux.f(nan_en);
            b_flux.p=deflux(:,:,nan_en);
            b_flux.f_units='eV';

            clear ind;
            if doPSD
                ind{1}=1:2;     % 0-30
                ind{2}=3:10;    % 30-150
                ind{2}=5:8;    % 60-120
                ind{3}=11:12;   % 150-180
            elseif doDEFlux
                
            else
            end
                

            for k=1:3
                flux{k}.t=b_flux.t;
                flux{k}.f=b_flux.f;
                flux{k}.p=double(squeeze(mean(b_flux.p(:,ind{k},:),2)));
                flux{k}.p_units=varunits;
                %flux{k}.f_units='keV'; % energy_units; % gives eV which is wrong
                flux{k}.f_units='eV'; % changed to eV to match with DEF
            end

            % Huishans recommendation
            aind=0;
            if 1 % 2*perp/(par+apar)
                aind=aind+1;
                an{aind}.p=log(2*flux{2}.p./(flux{1}.p+flux{3}.p));
                an{aind}.t=b_flux.t;
                an{aind}.f=b_flux.f; % saves in eV                
                an{aind}.p_units='log(2*per/\newline (par+par))';
                an{aind}.f_units='eV';%b_flux.f_units;
                an{aind}.varname=varname;
            end
            if 1 % par/apar
                aind=aind+1;
                an{aind}.p=log((flux{1}.p./flux{3}.p));
                an{aind}.t=b_flux.t;
                an{aind}.f=b_flux.f;                
                an{aind}.p_units='log(par/apar)';
                an{aind}.f_units='eV';%an{aind}.f_units=b_flux.f_units;
                an{aind}.varname=varname;
            end
    end
end

%% Check if there is sub count level data that should be put to NaN


%% Plot PSD data
nPanels = 3;
h = irf_plot(nPanels,'newfigure');
for oo = 1:3
    pcolor(h(oo),flux{oo}.t,flux{oo}.f,log10(flux{oo}.p)')
    shading(h(oo),'flat')
    hc(oo) = colorbar('peer',h(oo));
    caxis(h(oo),[-5 1])
    set(h(oo),'yscale','log')
    ylabel(h(oo),flux{oo}.f_units)
    ylabel(hc(oo),['log_{10} ' flux{oo}.p_units])
end
irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(h(nPanels),'usefig');

%% Plot anisotropy data
nPanels = 2;
h = irf_plot(nPanels,'newfigure');
irf_colormap('poynting')
for oo = 1:nPanels
    pcolor(h(oo),an{oo}.t,an{oo}.f,(an{oo}.p)')
    shading(h(oo),'flat')
    hc(oo) = colorbar('peer',h(oo));
    caxis(h(oo),2*[-1 1])
    set(h(oo),'yscale','log','ylim',[73 22000])
    ylabel(h(oo),an{oo}.f_units)
    ylabel(hc(oo),['log_{10} ' an{oo}.p_units])    
end
irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(h(nPanels),'usefig');
irf_pl_mark(h,toepoch([2007 08 31 10 17 30;2007 08 31 10 18 30])')
    