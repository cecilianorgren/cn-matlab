
varname = irf_ssub('C?_CP_PEA_PITCH_3DXH_PSD',sc);
[or_var,dobj,varmat,varunits] = c_caa_var_get(['Data__',varname]);
[or_bg,dobj_bg,varmat_bg,varunits_bg]=c_caa_var_get(['BackgroundLevel',xtra_vn,'__',varname]);

energy = varmat.dep_x{3}.data(1,:);
not_en_nan=~isnan(energy);
energy_units = c_caa_var_get(or_var.DEPEND_3,'units');
%%
var=or_var;
var_bg=or_bg;
var.data(or_var.data == or_var.FILLVAL) = NaN;
var.data(var_bg.data > var.data) = NaN;
%
% taking nanmean of azimuthal angles
% [time azimuth pitchangle energy] -> [time pitchangle energy]
data=squeeze(nanmean(var.data(:,:,:,:),2));

if 0
% first check to see if more than [certain vale] are NaNs, then we put that 
% value to NaN as well
data_nan = isnan(var.data);
data_nan_mean_azimuth = squeeze(mean(data_nan,2)); % [time pitchangle energy]
more_nan=find(data_nan_mean_azimuth>0.955);
data=squeeze(nanmean(var.data(:,:,:,:),2));
data(more_nan)=NaN;
end

data=data(:,:,not_en_nan);

% Definition of par, perp, apar
ind{1}=1:2;     % 0-30
ind{2}=5:8;    % 30-150
ind{3}=11:12;   % 150-180

for k=1:3
    flux{k}.t=varmat.t;
    flux{k}.f=energy(not_en_nan); %nan_en=~isnan(b_flux.f);
    flux{k}.p=double(squeeze(nanmean(data(:,ind{k},:),2)));
    flux{k}.p_units=varunits;    
    flux{k}.f_units='eV';
end

% Huishans recommendation
aind=0;
if 1
    aind=aind+1;
    an{aind}.p=log(flux{1}.p./(flux{3}.p));
    an{aind}.t=flux{aind}.t;
    an{aind}.f=flux{aind}.f; % saves in eV                
    an{aind}.p_units='log(par/apar)';
    an{aind}.f_units=flux{aind}.f_units;%b_flux.f_units;
    an{aind}.varname=varname;
    an{aind}.p(an{aind}.p>50)=NaN;
    an{aind}.p(an{aind}.p<-50)=NaN;
end
if 1
    aind=aind+1;
    an{aind}.p=log(2*flux{2}.p./(flux{1}.p+flux{3}.p));
    an{aind}.t=flux{aind}.t;
    an{aind}.f=flux{aind}.f; % saves in eV         
    an{aind}.p_units='log(2*per/\newline (par+par))';
    an{aind}.f_units=flux{aind}.f_units;%b_flux.f_units;
    an{aind}.varname=varname;
    an{aind}.p(isnan(an{1}.p))=NaN; % where the par apar is 'NaN' (too big/small)
    an{aind}.p(isnan(an{1}.p))=NaN; % -||-
end
