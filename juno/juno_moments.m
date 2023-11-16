% load('/Users/cecilia/Data/JUNO/Juno_data.mat')
addpath(genpath('/Users/cecilia/Data/JUNO/'))
units = irf_units;
mass = units.mp;
e = units.e;

it = 1:size(DATA_DIM1,1);
time = EpochTT(datestr(datetime(DIM0_UTC,'InputFormat','uuuu-DDD''T''HH:mm:ss.SSSSS'),'yyyy-mm-ddTHH:MM:SSZ'));
% DATA_DIM1
%            /* if background removed, VALID_MINIMUM, can be <0 */                   
%            VALID_MINIMUM     = -1.25e+14 /* = -999998/8e-09 */                     
%            VALID_MAXIMUM     =  2.81e+14 /* = 2250000/8e-09 */                     
%            MISSING_CONSTANT  = -999999                
%            MISSING_CONSTANT  = -999999                                             
%            UNIT              = "1/(m^2 sr s)"                                      
%            DESCRIPTION       = "DATA: Differential Energy Flux (SI units)          
%                                 64 Energy x 78 Look Directions.   
data_def = DATA_DIM1(it,:,:) + BACKGROUND_DIM1(it,:,:); 
% Is this above already the differential energy flux? Yes, seems so...


% Tried to see how doing these be below would affect the results
% data = -DATA_DIM1(it,:,:); % UNIT = "1/(m^2 sr s)"  ... Note: Flux, not psd. How to go from flux to psd?
% data(data<0) = 0;

% Baumjohann and Treuman, Basic space plasma physics:
% J = (v^2/m)*f
% Apply below


energy = DIM1_E_DIM1(it,:,:); % eV
energy_keV = energy/1e3; % eV
azim = DIM2_AZIMUTH_DESPUN_DIM1(it,:,:);
elev = DIM2_ELEVATION_DIM1(it,:,:); 
delta_elev = 11.25;
delta_azim = 22.50;
delta_energy_plus = energy;
diff_energy =  diff(energy,1,2);
energy_plus = energy; % will overwrite later, just to get dimensions right
energy_minus = energy;

energy_plus(:,1:end-1,:)  = energy(:,1:end-1,:) + 0.5*diff_energy;
energy_plus(:,end,:)      = energy(:,end,:) +     0.5*diff_energy(:,end,:);
energy_minus(:,1:end-1,:) = energy(:,1:end-1,:) - 0.5*diff_energy;
energy_minus(:,end,:)     = energy(:,end,:) -     0.5*diff_energy(:,end,:);
% plot(squeeze(azim(1,1,:)),squeeze(elev(1,1,:)),'*')

%data_def = data_flux/1e4./energy_keV; % 1/(m^2 s sr) -> 1/( cm^2 s sr keV), differential particle flux 
%data_dpf = data_flux/1e4./energy_keV; % 1/(m^2 s sr) -> 1/( cm^2 s sr keV), differential particle flux 
% iE = 1:64; semilogy(1:64,squeeze(energy),1:64,squeeze(energy_minus),1:64,squeeze(energy_plus))

data_def_cm = data_def*1e-4; 
data_def_tmp = data_def_cm;
data_def_tmp(isnan(data_def_tmp)) = 0;
data_def_omni = mean(data_def_tmp,3); % mean, but also need to divide 
ePDistOmni = PDist(time,data_def_omni,'omni',energy(:,:,1));
irf_spectrogram(ePDistOmni.specrec('energy'),'log'); hca = gca; hca.YScale = 'log';

%%
% McComas (JADE) Eq. 4.4: PSD = DEF/(2*E/m*q)^2
data_psd_m = data_def./(2*energy/mass/e).^(2); % phase density m^(-6)s^3 
data_psd_km = data_psd_m*1e18; % phase density km^(-6)s^3 
data_psd_omni_m = data_def_omni./(2*energy/mass/e).^(2); % phase density m^(-6)s^3 
data_psd_omni_km = data_psd_omni_m*1e18; % phase density km^(-6)s^3 

% data_psd_km = data_def_omni./energy.*0.53707.*((9.10956e-31./1.674e-27).^(2)); % phase density km^(-6)s^3 
% data_psd = data_psd_km*1e-18; % phase density m^(-6)s^3 


%delta_v_minus = sqrt(2*qe*energyupper/mass);
v_upper = sqrt(2*e*energy_plus/mass); % m/s
v_lower = sqrt(2*e*energy_minus/mass); % m/s
v_center = (v_upper + v_lower)/2;
v_delta = v_upper - v_lower;

solid_angle = sind(elev)*delta_elev*delta_azim*(pi/180)^2; % last bit transforms from angle to rad
v2dv = v_center.*v_delta.^2; % (m/s)^3
volume_element = v2dv.*solid_angle;


%data_psd = data_flux.*(mass./(v_center.^2)); % From flux to psd? s^3/m^6



dn = data_psd_m.*volume_element; % m^-3
n_psd_m3 = irf.nansum(irf.nansum(dn,3),2); % m^-3
n = n_psd_m3*1e-6; % cm^-3 NOT!, need to convert DATA from flux to psd, 
                       % maybe that will help, it was mentioned in the JUNO
                       % instruction book.
                       
                       
% For the oarticle flux we also need to add the direction vx, vy, vz
dvx = data_psd_m.*volume_element.*-cosd(elev).*sind(azim).*v_center;
dvy = data_psd_m.*volume_element.*-sind(elev).*sind(azim).*v_center;
dvz = data_psd_m.*volume_element.*-cosd(elev).*v_center;


vx = irf.nansum(irf.nansum(dvx,3),2)./n;
vy = irf.nansum(irf.nansum(dvy,3),2)./n;
vz = irf.nansum(irf.nansum(dvz,3),2)./n;


vx(abs(vx)>prctile(abs(vx),95))= NaN;
vy(abs(vy)>prctile(abs(vy),95))= NaN;
vz(abs(vz)>prctile(abs(vz),95))= NaN;

n(abs(n)>prctile(abs(n),95))= NaN;


tsN = irf.ts_scalar(time,n);
tsV = irf.ts_vec_xyz(time,[vx vy vz]);

%%
time_end = time.epochUnix;
Electron_caculate_data=[time_end,data_def_omni];%输入通量单位为1/(cm^2 s sr keV)
[Te,Ne,ePSD,Ve]=Juno_neTe_spart_keV(energy_keV(:,:,1),Electron_caculate_data);
specrec.p = ePSD(:,2:end);
specrec.t = ePSD(:,1);
specrec.f = energy_eV;
specrec.p_label = {''};
%%
h = irf_plot(7);

hca = irf_panel('n 123');
irf_plot(hca,tsN)

hca = irf_panel('v 123');
irf_plot(hca,tsV)

hca = irf_panel('omni 0');
irf_spectrogram(hca,specrec)
hca.YScale = 'log';
%colorbar(hca)

hca = irf_panel('omni 123');
pdist = ePDistOmni;
irf_spectrogram(hca,pdist.specrec('energy'))
hca.YScale = 'log';

hca = irf_panel('omni 1');
pdist = ePDistOmni(1:3:ePDistOmni.length);
irf_spectrogram(hca,pdist.specrec('energy'))
hca.YScale = 'log';

hca = irf_panel('omni 2');
pdist = ePDistOmni(2:3:ePDistOmni.length);
irf_spectrogram(hca,pdist.specrec('energy'))
hca.YScale = 'log';

hca = irf_panel('omni 3');
pdist = ePDistOmni(3:3:ePDistOmni.length);
irf_spectrogram(hca,pdist.specrec('energy'))
hca.YScale = 'log';


%irf_plot({tsN.resample(time(1:3:time.length)),tsVx.resample(time(1:3:time.length))})
hlinks = linkprop(h,{'XLim'});
irf_plot_axis_align(h)
%% Notes
%      OBJECT              = CONTAINER                                               
%        NAME              = DATA_DIM1                                               
%        START_BYTE        = 321                                                     
%        BYTES             = 312 /* = 78 * 4-bytes */                                
%        REPETITIONS       = 64                                                      
%        DESCRIPTION       = "DATA_DIM1,                                             
%                             2D array of data, 1st and 2nd Dimensions."             
%                                                                                    
%        OBJECT              = CONTAINER                                             
%          NAME              = DATA_DIM2                                             
%          START_BYTE        = 1                                                     
%          BYTES             = 4                                                     
%          REPETITIONS       = 78                                                    
%          DESCRIPTION       = "DATA_DIM2,                                           
%                               1D array of data, 2nd Dimension."                    
%                                                                                    
%          OBJECT              = COLUMN                                              
%            NAME              = DATA                                                
%            DATA_TYPE         = PC_REAL                                             
%            START_BYTE        = 1                                                   
%            ITEMS             = 1                                                   
%            ITEM_BYTES        = 4                                                   
%            BYTES             = 4                                                   
%            /* if background removed, VALID_MINIMUM, can be <0 */                   
%            VALID_MINIMUM     = -1.25e+14 /* = -999998/8e-09 */                     
%            VALID_MAXIMUM     =  2.81e+14 /* = 2250000/8e-09 */                     
%            MISSING_CONSTANT  = -999999                                             
%            UNIT              = "1/(m^2 sr s)"                                      
%            DESCRIPTION       = "DATA: Differential Energy Flux (SI units)          
%                                 64 Energy x 78 Look Directions.                    
%                                 "                                                  
%      /* RJW, DATA, f, 2, 64, 78 */                                                 
%          END_OBJECT          = COLUMN                                              
%        END_OBJECT          = CONTAINER                                             
%      END_OBJECT          = CONTAINER                                               
%                                                                                    
%      OBJECT              = CONTAINER                                               
%        NAME              = DATA_SIGMA_DIM1                                         
%        START_BYTE        = 20289                                                   
%        BYTES             = 312 /* = 78 * 4-bytes */                                
%        REPETITIONS       = 64                                                      
%        DESCRIPTION       = "DATA_SIGMA_DIM1,                                       
%                             2D array of data, 1st and 2nd Dimensions."             
%                                                                                    
%        OBJECT              = CONTAINER                                             
%          NAME              = DATA_SIGMA_DIM2                                       
%          START_BYTE        = 1                                                     
%          BYTES             = 4                                                     
%          REPETITIONS       = 78                                                    
%          DESCRIPTION       = "DATA_SIGMA_DIM2,                                     
%                               1D array of data, 2nd Dimension."                    
%                                                                                    
%          OBJECT              = COLUMN                                              
%            NAME              = DATA_SIGMA                                          
%            DATA_TYPE         = PC_REAL                                             
%            START_BYTE        = 1                                                   
%            ITEMS             = 1                                                   
%            ITEM_BYTES        = 4                                                   
%            BYTES             = 4                                                   
%            VALID_MINIMUM     =  0                                                  
%            VALID_MAXIMUM     =  1.25e+13 /* = 100000/8e-09 */                      
%            MISSING_CONSTANT  = -999999                                             
%            UNIT              = "1/(m^2 sr s)"                                      
%            DESCRIPTION       = "DATA_SIGMA                                         
%                                 1-sigma uncertainties on values in object DATA,    
%                                 such that true value = DATA +/- DATA_SIGMA.        
%                                 See DATA entry above for size information."    
