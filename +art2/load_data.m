cd /Users/Cecilia/Data/BM/20070831
t1=[2007 08 31 10 13 00];
t2=[2007 08 31 10 25 00];
tint=toepoch([t1;t2])';

% Just load everything, it's easier
fig1=1; fig2=1; fig3=1; fig4=1;

%% Magnetic field
if 1;%any([fig1 fig2])
    c_eval('caa_load(''C?_CP_FGM_FULL'');',1:4);
    c_eval('gseB?=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');',1:4);
    c_eval('gsmB?=c_coord_trans(''gse'',''gsm'',gseB?,''cl_id'',?);',1:4);
    c_eval('diB?=c_coord_trans(''gse'',''isr2'',gseB?,''cl_id'',?);',1:4);
    if ~exist('diB3','var'); load diB; end
end
if 1;%1, % read STAFF data form all sc
    load mBS
    tint=[dBS3(1,1) dBS3(end,1)];
    c_eval('B?staff=c_coord_trans(''DSC'',''gsm'',dBS?,''cl_id'',?);',3:4);
    c_eval('distaB?=c_coord_trans(''DSC'',''isr2'',dBS?,''cl_id'',?);',3:4);
    c_eval('gsmBfs?=c_fgm_staff_combine(gsmB?(:,1:4),B?staff);',3:4)
    c_eval('gseB?staff=c_coord_trans(''DSC'',''gse'',dBS?,''cl_id'',?);',3:4);
    c_eval('gseBfs?=c_coord_trans(''gsm'',''gse'',gsmBfs?,''cl_id'',?);',3:4);
    c_eval('diBfs?=c_coord_trans(''gsm'',''isr2'',gsmBfs?,''cl_id'',?);',3:4);
end
%% Ion velocity
if any([fig1])
    sclist = 3;
    c_eval('caa_load(''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',sclist);
    c_eval('gseVhia?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sclist);
    c_eval('gsmVhia?=c_coord_trans(''gse'',''gsm'',gseVhia?,''cl_id'',?);',sclist);    
    c_eval('caa_load(''C?_CP_CIS-CODIF_HS_He1_MOMENTS'');',sclist);
    c_eval('gseVcod?=getmat(C?_CP_CIS_CODIF_HS_He1_MOMENTS,''velocity__C?_CP_CIS_CODIF_HS_He1_MOMENTS'');',sclist);
    c_eval('gsmVcod?=c_coord_trans(''gse'',''gsm'',gseVcod?,''cl_id'',?);',sclist);    
    
end


%% Spacecraft location
if fig2
    c_eval('gsePos?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);         
    %c_eval('diPos?=c_coord_trans(''gse'',''isr2'',gsePos?,''cl_id'',?);',3:4);         
    %c_eval('gsmPos?=c_coord_trans(''gse'',''gsm'',cn_toepoch(t1,t2,gsePos?),''cl_id'',?);',3:4)
end
%% Spacecraft potential
if fig1
    c_eval('scP?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',3:4);         
    
end


%% Coordinate system
% It is better to have a common coordinate system for the two satellites,
% therefore I should use gse-coordinates instead of isr2 (which are
% individual for each satellite, since they are slightly tilted). By uing
% this system, it becomes impossible to have the perpendicular electric
% field exctly in the spin plane, but the difference should be small.
% The angle between the two spin-axes were 2 deg, so any differences can be
% neglected.
if any([fig1 fig2 fig3]) % spin plane etc
    bfield = irf_add(0.5,irf_norm(gseB3),0.5,irf_norm(gseB4)); % along B-field
    c_eval('diz? = [diB?(:,1) repmat([0 0 1],size(diB?,1),1)];',3:4);
    c_eval('gsez? = c_coord_trans(''isr2'',''gse'',diz?,''cl_id'',?);',3:4);
    spinaxis = irf_norm(irf_add(0.5,gsez3,0.5,gsez4));
    spinplane = irf_cross(spinaxis,bfield); % in spin plane
    
    % time series units vectors in gse-coordinates
    third = bfield; % along magnetic field
    first = spinplane; % in spin plane, perpendicular to magnetic field
    second = irf_cross(third,first);
    
end

if any([fig2 fig3]) % spacecraft location
    d0 = irf_add(0.5,gsePos3,0.5,gsePos4); % vector to center between
    dpos = irf_add(1,gsePos3,-1,gsePos4); % vector from C4 to C3
    dposper2 = irf_norm(irf_cross(bfield,dpos));
    dposper1 = irf_norm(irf_cross(bfield,dposper2));
    position = dposper1;
end

if 1 % difference between spin plane and R_C3 - R_C4
    dotprod = irf_dot(third,position);
end

%% Spacecraft location
if any([fig2 fig3])
    %c_eval('par? = irf_norm(gseB?);',3:4);
    %par = irf_add(0.5,par3,0.5,par4);
    % Get sc positions in gse-system (not ISR2 because it moves and is 
    % different between c3 and c4).
    dpos = irf_add(1,gsePos3,-1,gsePos4);
    dp1 = irf_dot(first,dpos);  % in spin plane, perp to B  
    dp2 = irf_dot(second,dpos); % BxSP
    dp3 = irf_dot(third,dpos);  % par to B
    dpar = dp3; % since first is along b-field
    dper = irf_add(1,dpos,-1,dpar);
    
    %didper = c_coord_trans('gse','isr2',dper,'cl_id',)
    %dz=irf_dot(irf_add(1,gsePos4,-1,gsePos3),z);
    %dzav=sum(dz(:,2)/size(dz,1));
    %dy=irf_dot(irf_add(1,gsePos4,-1,gsePos3),y); % y4-y3
    %dx=irf_dot(irf_add(1,gsePos4,-1,gsePos3),x); % x4-x3
    %dyav=mean(dy(:,2));
    %dxav=mean(dx(:,2));
    
    %dpar=[dxav dyav dzav]*x(2:4)';
    %facdy=[dxav dyav dzav]*y(2:4)';
    %facdz=[dxav dyav dzav]*z(2:4)';
end

%% Electric field
if any([fig1 fig2])
    if ~exist('diE3','var'); 
        load diE; diE3(:,4)=0; diE4(:,4)=0; 
        c_eval('gseE? = c_coord_trans(''isr2'',''gse'',diE?,''cl_id'',?);',3:4);
    end
    if ~exist('diB3','var'); load diB; end

    lim_ang = 45;
    
    % construct parallel electric field
    c_eval('[diEtemp? eb_anglong?]=irf_edb(irf_tlim(diE?,tint+[-1 1]*60),irf_tlim(diB?(:,1:4),tint+[-1 1]*60),lim_ang,''Epar'');',3:4);
    c_eval('eb_anglong? = [diEtemp?(:,1) eb_anglong?];',3:4);
    %c_eval('diz? = [diB?(:,1) repmat([0 0 1],size(diB?,1),1)];',3:4);
    %c_eval('dib? = irf_norm(diB?);',3:4)
    %c_eval('diper? = irf_cross(diz?,dib?);',3:4)
    %c_eval('dipar? = irf_cross(diper?,diz?);',3:4)
    c_eval('Epar?=irf_dot(irf_tlim(diE?,tint),third);',3:4); % along B
    c_eval('Eper?=irf_dot(irf_tlim(diE?,tint),first);',3:4); % in spin plane perp to Bs
    % the last direction is not interesting since its just the other
    % projection of the field that is along B.

    fmin = 0.5; fmax = 0; order = 5;
    c_eval('EparAC? = irf_filt(Epar?,fmin,fmax,cn.f(Epar?),order);',3:4);
    c_eval('EparDC? = irf_add(1,Epar?,-1,EparAC?);',3:4);
end

%% Wave magnetic field
    c_eval('Bpar?=irf_dot(irf_tlim(diBfs?,tint),third);',3:4); % along B
    c_eval('Bper?=irf_dot(irf_tlim(diBfs?,tint),first);',3:4); % in spin plane perp to Bs
    c_eval('staBpar?=irf_dot(irf_tlim(distaB?,tint),third);',3:4); % along B
    c_eval('staBper1?=irf_dot(irf_tlim(distaB?,tint),first);',3:4); % in spin plane perp to Bs
    c_eval('staBper2?=irf_dot(irf_tlim(distaB?,tint),second);',3:4); % in spin plane perp to Bs
    % the last direction is not interesting since its just the other
    % projection of the field that is along B.

    fmin = 5; fmax = 0; order = 5;
    c_eval('BparAC? = irf_filt(Bpar?,fmin,fmax,cn.f(Bpar?),order);',3:4);
    c_eval('BparDC? = irf_add(1,Bpar?,-1,BparAC?);',3:4);

%% Electron temperature
c_eval('parTe?=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
%c_eval('peaTepar?z=cn_toepoch(t1,t2,parTe?);',3:4);
%    Teparav=mean([peaTepar3z(:,2);peaTepar4z(:,2)])*8.61734*10;
%    c_eval('peaTeper?z=cn_toepoch(t1,t2,perTe?);',3:4);
%    Teperav=mean([peaTeper3z(:,2);peaTeper4z(:,2)])*8.61734*10;
%    Teav=(Teparav+Teperav)/2;

    

%% Anisotropy data   
if fig1
    loadfrommat=1;
    doDEFlux = 0;
    doDPFlux = 0;
    doPSD = 1;
    
    switch loadfrommat
        case 1
            load an3 an4
        case 0
            %%
            sc=4;
            art2.do_anisotropy;
            save('an4','an')
            sc=3;
            art2.do_anisotropy;
            save('an3','an')
            if 0
            sc = 4;
            if doPSD; varname=irf_ssub('C?_CP_PEA_PITCH_3DXH_PSD',sc);
            elseif doDEFlux; varname=irf_ssub('C?_CP_PEA_3DXPH_DEFlux',sc);
            else varname=irf_ssub('C?_CP_PEA_3DXPH_DPFlux',sc); end

            xtra_vn=''; % could be _HEEA or _LEEA when using pitch full
            [or_var,dobj,varmat,varunits]=c_caa_var_get(['Data',xtra_vn,'__',varname]);
            varname2=irf_ssub(['Sweep_Energy',xtra_vn,'__',varname],sc);
            %varname3=irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
            energies=c_caa_var_get(varname2,'mat');
            energy_units=c_caa_var_get(varname2,'units'); % These give eV, should be keV
            %piichangle=c_caa_var_get(varname3,'mat');

            % Taking away all fluxes that are below one, these are something else,
            % what?
            var=or_var;
            var.data(or_var.data == or_var.FILLVAL)=NaN;

            % taking mean of azimuthal angles
            % [time azimuth pitchangle energy]
            deflux=squeeze(nanmean(var.data(:,:,:,:),2));
            % [time pitchangle energy]

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
                ind{3}=11:12;   % 150-180
            elseif doDEFlux
                
            else
            end
                

            for k=1:3
                flux{k}.t=b_flux.t;
                flux{k}.f=b_flux.f; % converts to eV
                flux{k}.p=double(squeeze(nanmean(b_flux.p(:,ind{k},:),2)));
                flux{k}.p_units=varunits;
                flux{k}.f_units='keV'; % energy_units; % gives eV which is wrong
                flux{k}.f_units='eV'; % changed to eV to match with DEF
            end

            % Huishans recommendation
            aind=0;
            if 1
                aind=aind+1;
                an{aind}.p=log(flux{2}.p./(flux{3}.p));
                an{aind}.t=b_flux.t;
                an{aind}.f=b_flux.f; % saves in eV                
                an{aind}.p_units='log(2*per/\newline (par+par))';
                an{aind}.f_units='eV';%b_flux.f_units;
                an{aind}.varname=varname;
            end
            if 1
                aind=aind+1;
                an{aind}.p=log(2*flux{2}.p./(flux{1}.p+flux{3}.p));
                an{aind}.t=b_flux.t;
                an{aind}.f=b_flux.f;                
                an{aind}.p_units='log(par/apar)';
                an{aind}.f_units='eV';%an{aind}.f_units=b_flux.f_units;
                an{aind}.varname=varname;
            end
            end
    end
end

%% Density
if fig1
    if ~exist('peaNe3','var'); load matlabN; end
end