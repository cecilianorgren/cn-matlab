function results = cn__pitch_full_subspin(sc)
xtra={'_HEEA','_LEEA'};
for sensor=1:2;    
    varname=irf_ssub('C?_CP_PEA_PITCH_FULL_DEFlux',sc);
    xtra_vn=xtra{sensor}; % could be _HEEA or _LEEA when using pitch full
    [variable,dataobject,varmat,dataunits]=c_caa_var_get(['Data',xtra_vn,'__',varname]);
    %varname2=irf_ssub(['Sweep_Energy',xtra_vn,'__',varname],sc);
    %varname3=irf_ssub(['Sweep_Azimuth',xtra_vn,'__',varname],sc);
    %varname3=irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
    %energies=c_caa_var_get(varname2,'mat');
    %energy_units=c_caa_var_get(varname2,'units'); % These give eV, should be keV
    
    phivar=c_caa_var_get(variable.DEPEND_1);phivar=fillval_to_nan(phivar);phi=phivar.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
    envar=c_caa_var_get(variable.DEPEND_3);envar=fillval_to_nan(envar);en=envar.data(1,:);nan_en=isnan(en);en(nan_en)=[];
    pitchvar=c_caa_var_get(variable.DEPEND_2); pitch_angle=pitchvar.data(1,1:size(pitchvar.data,2))
    enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
    enlabel=enunits;
    variable.data(:,:,:,nan_en)=[]; % remove NaN energy data
    variable.data(:,nan_phi,:,:)=[]; % remove NaN energy data
    %variable.data=permute(variable.data,[1 3 2 4]); % permute in order [time, polar, azimuth, energery]
    pitchangle=variable.data; 
    t=varmat.t(:);
    dataraw=variable.data;
    dataraw(dataraw == variable.FILLVAL)=NaN; % FILLVALs put to NaN
    [en,ix]=sort(en); % sort energy in ascending order
    dataraw=dataraw(:,:,:,ix); % sort data accordingly
    dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
    data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
    tt=repmat(t(:),1,length(phi));
    phiphi=tt;
    [tt]=subspintime(dataobject,phi);
    variable=fillval_to_nan(variable); % FILLVALs put to NaN
    ind_data = data; ind_data(~isnan(data)) = 1; ind_data(isnan(data)) = 0;
data_with_nan=data;
    
   res.tt = tt;                    % time axis
%res.tt_deltaplus = [];
%res.tt_deltaminus = [];
res.en = en;                    % energy levels
res.phi = phi;                  % azimuthal angles
%res.theta = theta;              % pitch angles
res.data=data_with_nan;
%res.dtsampling=dtsampling;
%res.omni=data_omni;
res.pitch_angle=pitch_angle;
%res.phiphi=phiphi;
res.enlabel=enlabel;
res.dataunits=dataunits;
results{sensor}=res;
end
end

%----------------------------------------------------------
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
spin_period=double(dtplus+dtminus);
dtsampling=spin_period/length(phi);
for j=length(phi):-1:1,
    tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
end
tt=reshape(tt',numel(tt),1);
end

%----------------------------------------------------------
function out=fillval_to_nan(in)
% FILLVALs put to NaN
out=in;
out.data(in.data == in.FILLVAL)=NaN;
end     