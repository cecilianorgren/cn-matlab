function specrec = c_caa_distribution_data(product,detector)
% C_CAA_DISTRIBUTION_DATA Prepares data to be plotted with
% c_caa_distribution_function.
%   data_structure = C_CAA_DISTRIBUTION_DATA(data_product,detector)
%       data_product - for example: 'C3_CP_PEA_3DXPH_PSD'
%       detector - 'LEEA' or 'HEEA', needs to be specified when using data
%                  product ' _FULL'
%   
%
%    e_psd3 = 
% 
%         product: 'C3_CP_PEA_3DXPH_PSD'
%               t: [69344x1 double]
%              dt: 0.0646
%        en_label: 'Energy [eV]'
%           en_cs: [26x1 single]
%          en_pol: [27x1 single]
%         f_label: 'Pitch angle'
%            f_cs: [15 45 75 105 135 165]
%           f_pol: [0 30 60 90 120 150 180]
%         p_label: 'Log PSD [s^3/km^6]'
%               p: [69344x6x26 double]
%            p_bg: [69344x6x26 double]
%
%   Examples:
%       data_structure = C_CAA_DISTRIBUTION_DATA('C3_CP_PEA_3DXPH_PSD');
%       h=c_caa_distribution_function(tint,'polar',data_structure);
%
% See also cn_plot_distribution_function_combined.m,
%          testing_distribution_function_combined.m
distr=product(max(strfind(product,'_'))+1:end);
specrec.product=product;

if any(strfind(product,'PITCH')) % Something fishy with background level
    disp('Something fishy with background level...')
    if 0;any(strfind(product,'FULL')) % contains both detectors
        detector = input('Which detector? LEEA (1) or HEEA (2)?');        
        switch detector
            case 1; detector_str='_LEEA';
            case 2; detector_str='_HEEA';
        end
        specrec.detector=detector_str(2:end);
    else
        detector_str='';
    end
    
    if any(strfind(product,'FULL')); detector_str=['_',detector];
    else detector_str=''; 
    end
    
    % Load pitch angle data
    [caaData,dataobject,data,data_units]=c_caa_var_get(['Data',detector_str,'__',product]);
    [caabg,~,bg,bg_units]=c_caa_var_get(['BackgroundLevel',detector_str,'__',product]);    

    % Pitch angles
    theta=data.dep_x{2}.data(1,:);
    theta_plus=data.dep_x{2}.df.plus;
    theta_minus=data.dep_x{2}.df.minus;    
    nan_theta=isnan(theta);theta(nan_theta)=[];
    theta_pol=[theta-theta_minus theta(end)+theta_plus(end)]; 
    theta_cs=theta;
    
    % Energy levels
    en=data.dep_x{3}.data(1,:);
    en_plus=data.dep_x{3}.df.plus;
    en_minus=data.dep_x{3}.df.minus;    
    nan_en=isnan(en);en(nan_en)=[];en_plus(nan_en)=[];en_minus(nan_en)=[];                    
    en_pol=[en+en_plus en(end)-en_minus(end)];
    en_cs=en;

    % Sub spin angles
    phi=data.dep_x{1}.data(1,:);
    nan_phi=isnan(phi);phi(nan_phi)=[];

    % Construct subspin resolution data
    dataraw=data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
    dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
    data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
    [tt dtsampling]=subspintime(dataobject,phi); % define subspin time from angle phi.

    % Same for background level
    bgraw=bg.data; bgraw(:,:,:,nan_en)=[];bgraw(:,nan_phi,:,:)=[];
    bgraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
    bg=reshape(bgraw,size(bgraw,1)*size(bgraw,2),size(bgraw,3),size(bgraw,4));

    % Store all data in specrec    
    specrec.t=tt; % time
    specrec.dt=dtsampling/2;        
    
    specrec.en_label='Energy [eV]'; % energy
    [specrec.en_cs en_cs_order]=sort(en_cs);
    [specrec.en_pol en_pol_order]=sort(en_pol);
       
    specrec.f_label='Pitch angle'; % pitch angle
    specrec.f_cs=theta_cs;
    specrec.f_pol=theta_pol;
    
    specrec.p_label=['Log ' distr ' [' data_units ']']; % data
    specrec.p=double(data(:,:,en_cs_order)); % average over time
    specrec.p_bg=double(bg(:,:,en_cs_order)); % average over time        
end
if any([strfind(product,'3DXPH')])
    % Load pitch angle data
    res=c_caa_construct_subspin_res_data(['Data__', product]);
    bg=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);
    
    SEDL=c_caa_var_get(['Sweep_Energy_DeltaLower__', product],'mat');
    SEDU=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product],'mat');
    en_nan=isnan(SEDL(1,:));SEDL(:,en_nan)=[]; [SEDL SEDL_order]=sort(SEDL(1,2:end));
    en_nan=isnan(SEDU(1,:));SEDU(:,en_nan)=[]; [SEDU SEDU_order]=sort(SEDU(1,2:end)); 
    en_pol=[res.en-SEDL res.en(end)+SEDU(end)];    
    
    dtheta=diff(res.theta(1:2))/2;
    theta_pol=[res.theta-dtheta res.theta(end)+dtheta];
    
    % Store all data in specrec    
    specrec.t=res.tt;
    specrec.dt=res.dtsampling/2;
    
    specrec.en_label='Energy [eV]'; % energy
    specrec.en_cs=res.en(:);
    specrec.en_pol=en_pol(:);
    
    specrec.f_label='Pitch angle'; % pitch angle
    specrec.f_cs=res.theta;
    specrec.f_pol=theta_pol;
    
    %specrec.p=res.pitch_angle(:,:);
    specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   
    specrec.p=res.data(:,:,:);
    specrec.p(find(specrec.p==0))=NaN; 
    specrec.p_bg=bg.data(:,:,:);        
end
if any([strfind(product,'CODIF_HS'),...
        strfind(product,'CODIF_LS'),...
        strfind(product,'HIA_LS'),...
        strfind(product,'HIA_HS')]);
    
    % Load pitch angle data
%     in case this doesnt work, change to something like:
%     caaSEDL=c_caa_var_get(['delta_plus_energy_table__', products{k}]);
%     caaSEDU=c_caa_var_get(['delta_minus_energy_table__', products{k}]);                            
%     SEDL=flipdim(caaSEDL.data(1,:),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
%     SEDU=flipdim(caaSEDU.data(1,:),2); en_nan=isnan(SEDU);SEDU(en_nan)=[];
    res=c_caa_construct_subspin_res_data(['x3d_ions__', product]);
    SEDU=c_caa_var_get(['delta_plus_energy_table__', product],'mat'); 
    SEDL=c_caa_var_get(['delta_minus_energy_table__', product],'mat'); 
    SEDU=torow(sort(SEDU(1:31))); SEDL=torow(sort(SEDL(1:31)));
    en_pol=[res.en-SEDL res.en(end)+SEDU(end)];
    
    dtheta=diff(res.theta(1:2))/2;
    theta_pol=[res.theta-dtheta res.theta(end)+dtheta];
    
    % No background level data?

    % Store all data in specrec  
    specrec.t=res.tt;
    specrec.dt=res.dtsampling/2;
    
    specrec.en_label='Energy [eV]'; % energy
    specrec.en_cs=res.en(:);
    specrec.en_pol=en_pol(:);
    
    specrec.f_label='Pitch angle'; % pitch angle
    specrec.f_cs=res.theta;
    specrec.f_pol=theta_pol;
        
    specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   
    specrec.p=res.data; 
    specrec.p_bg=[]; 
end

if 0 % still to fix
if any(strfind(product,'PAD'))
    [caaData,dataobject,Data,Data_units]=c_caa_var_get(['PAD_Electron_Dif_flux__',product]);

    theta=Data.dep_x{2}.data(1,:);
    nan_theta=isnan(theta);theta(nan_theta)=[]; % for data

    en=Data.dep_x{1}.data(1,:);
    nan_en=isnan(en);en(nan_en)=[]; % for data

    % Construct subspin reoslution data
    dataraw=Data.data; dataraw(:,nan_en,:)=[];

    % Pick out defined time interval
    dtsampling=(Data.t(2)-Data.t(1))/2;
    tmin=Data.t-dtsampling;
    tmax=Data.t+dtsampling;
    tind=find(tmin<tint(2));
    tind=find(tmax(tind)>tint(1));

    specrec.f=theta;
    specrec.en=en'; % stigande
    specrec.t=Data.t(tind);
    specrec.p_label=['Log ' distr ' [' Data_units ']'];
    specrec.f_label='Pitch angle';
    specrec.data=squeeze(nanmean(permute(double(dataraw(tind,:,:)),[1 3 2]))); % order: time, pitch angle, energy
    data_0=find(specrec.data==0);
    specrec.data(data_0)=NaN;  
    specrec.bg=specrec.data;
    specrec.bg(:)=NaN;
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
end