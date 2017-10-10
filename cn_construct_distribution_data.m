function specrec = cn_construct_distribution_data(product,plot_type)
%
%
% See also cn_plot_distribution_function
distr=product(max(strfind(product,'_'))+1:end);
specrec.product=product;

switch plot_type
    case 'cross-section'
        if any(strfind(product,'PITCH')) % Something fishy with background level
            disp('Something fishy with background level...')
            % Load pitch angle data
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',product]);
            [caabg,~,bg,bg_units]=c_caa_var_get(['BackgroundLevel__',product]);

            % Pitch angles
            theta=Data.dep_x{2}.data(1,:);
            nan_theta=isnan(theta);theta(nan_theta)=[];
            % Energy levels
            en=Data.dep_x{3}.data(1,:);
            nan_en=isnan(en);en(nan_en)=[];
            % Sub spin angles
            phi=Data.dep_x{1}.data(1,:);
            nan_phi=isnan(phi);phi(nan_phi)=[];

            % Construct subspin resolution data
            dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            tt=subspintime(dataobject,phi); % Define subspin time from angle phi.

            % Same for background level
            bgraw=bg.data; bgraw(:,:,:,nan_en)=[];bgraw(:,nan_phi,:,:)=[];
            bgraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            bg=reshape(bgraw,size(bgraw,1)*size(bgraw,2),size(bgraw,3),size(bgraw,4));

            specrec.f=theta;
            [specrec.en en_order]=sort(en);
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' Data_units ']'];
            specrec.f_label='Pitch angle';
            specrec.en_label='eV';
            specrec.bg=double(bg(:,:,en_order)); % average over time
            specrec.data=double(data(:,:,en_order)); % average over time
        end
        if any([strfind(product,'3DXPH')])
            % Load pitch angle data
            res=c_caa_construct_subspin_res_data(['Data__', product]);
            bg=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);

            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(:,:);
            specrec.en=res.en(:);
            specrec.data=res.data(:,:,:);
            specrec.bg=bg.data(:,:,:);
            specrec.bg_en=bg.en(:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   

            %specrec.data=squeeze(nanmean(specrec.data,1));
            %specrec.bg=squeeze(nanmean(specrec.bg,1));
            energy=specrec.en;
            specrec.type=plot_type;
        end
        if any([strfind(product,'CODIF_HS'),...
                strfind(product,'CODIF_LS'),...
                strfind(product,'HIA_LS'),...
                strfind(product,'HIA_HS')])
            % Load pitch angle data
            res=c_caa_construct_subspin_res_data(['x3d_ions__', product]);

            % No background level data?

            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(:,:);
            specrec.en=res.en(:);
            specrec.data=res.data(:,:,:);
           % specrec.bg_en=bg.en(:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   

            %specrec.data=squeeze(nanmean(specrec.data,1));
            %specrec.bg=squeeze(nanmean(specrec.bg,1));
            %specrec.bg=specrec.data;
            %specrec.bg(:)=NaN;
            energy=specrec.en;
        end
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
    case 'polar'
       if any(strfind(product,'PITCH'))
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data_LEEA__',product]);

            % Set up grid
            theta_plus=Data.dep_x{2}.df.plus;
            theta_minus=Data.dep_x{2}.df.minus;
            theta=Data.dep_x{2}.data(1,:);
            nan_theta=isnan(theta);theta(nan_theta)=[]; % for data
            theta=[theta-theta_minus theta(end)+theta_plus(end)]; % grid

            en_plus=Data.dep_x{3}.df.plus;
            en_minus=Data.dep_x{3}.df.minus;
            en=Data.dep_x{3}.data(1,:);
            nan_en=isnan(en);en(nan_en)=[]; % for data
            en_plus(nan_en)=[];en_minus(nan_en)=[];                    
            en=[en+en_plus en(end)-en_minus(end)];

            phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];

            % Construct subspin reoslution data
            dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));

            [tt dtsampling]=subspintime(dataobject,phi); % Define subspin time from angle phi.


            specrec.p=theta;
            specrec.f=theta;
            specrec.en=flipdim(en,2); % stigande
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' Data_units ']'];
            specrec.f_label='Pitch angle';
            specrec.data=double(data(:,:,:));
            specrec.data=flipdim(specrec.data,3);
            data_0=find(specrec.data==0);
            specrec.data(data_0)=NaN; 
            specrec.en_label='Energy [eV]';
            specrec.type=plot_type;
            
       end
       if any(strfind(product,'3DXPH'))
                res=c_caa_construct_subspin_res_data(['Data__', product]);
                [caaSEDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', product]);
                [caaSEDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product]);
                SEDL=flipdim(SEDL(1,2:end),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
                SEDU=flipdim(SEDU(1,2:end),2); en_nan=isnan(SEDU);SEDU(en_nan)=[];   
    
                %[~,ind]=irf_tlim(res.tt,tint);
                
                specrec.t=res.tt;
                specrec.dt=res.dtsampling/2;
                dtheta=(res.theta(2)-res.theta(1))/2;
                specrec.f=[res.theta-dtheta res.theta(end)+dtheta];
                specrec.f_label='Pitch angle';
                specrec.p=res.pitch_angle(:,:);
                specrec.en=res.en(:);
                specrec.en=[res.en-SEDL res.en(end)+SEDU(end)];
                specrec.en_label='Energy [eV]';
                specrec.data=res.data(:,:,:);
                specrec.p_label={['Log ' distr ' [' res.dataunits ']']};

                data=find(specrec.data==0);
                specrec.data(data)=NaN; 
                specrec.type=plot_type;
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