function [energy,PA0,PA90,PA180]=c_caa_distribution_function(varargin)
% C_CAA_DISTRIBUTION_FUNCTION  Plot particle distribution.
%
% Plots a polar, cross-section or ... distribution of particle data.
% Time interval is specified as 'tint',[tint_in_epoch].
% Axis handle can be specified.
% Can plot PEACE PITCH and 3DXPH products.
% The type of plot is specified as an option: 
%   'polar' - Default, can plot one or two products at a time.
%   'cross-section' - Plots cross-sections of the distribution 
%                     function at 0, 90 and 180 degrees as well
%                     as zero-count level. Can plot one product.
%                     at a time.
% 
% 
% Example: c_caa_distributon_function(h,'C4_CP_PEA_PITCH_3DXH_DEFlux','tint',tint)
%          c_caa_distributon_function('C4_CP_PEA_3DXPH_PSD','tint',tint)
%          c_caa_distributon_function(h,'C4_CP_PEA_PITCH_3DXH_DPFlux','tint',tint,'cross-section')

% Current issues: o 3DRL: problem creating subspin time since sampling
% period seems very irregular. dt varies between 20 and 12 whereas it is
% constantly ~4 for 3DXH. 
% o Colorbar: While it is possible to plot two different products in the
% same plot, the color coding can not be separated. Don't know how to solve
% this.

% For PEACE only PITCH and 3DXPH because 3DXPH is the only one included 
% in c_caa_construct_subspin_res_data. The angular resolution also depends
% on this function.

% Matlab function surf is used. The grid is defined by the energy and pitch
% angle intervals, i.e. one value extra per dimension. The color of each
% square/interval is taken from the data.

% General figure options
set(0,'defaultAxesFontSize',14);
set(0,'defaultTextFontSize',14);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');

% Read input
[ax,args,nargs] = axescheck(varargin{:});
original_args=args;

% Default values that one can override with options.
plot_type='polar';
n_products=0;
emin=[]; % Defines the origin
products={};

% Check for products
l=1;
while l<=nargs
    if any(strfind(args{l},'RAP'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='RAP';
    elseif any(strfind(args{l},'CIS'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='CIS'; 
    elseif any(strfind(args{l},'PEA'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='PEA';       
    else
        switch lower(args{l})
            case 'tint'
                l=l+1;
                tint=args{l};
            case 'polar'
                plot_type='polar';
            case 'cross-section'
                plot_type='cross-section';              
            case 'emin'
                l=l+1;
                emin=args{l};
        end
    end
l=l+1;
end

% Return if not enough input is given.
% Reduce products if too much input is given.
if isempty(tint)
    disp('No time interval was given!')
    return;
end
switch n_products
    case 0
        disp('No particle product was given!')
        return;
    case 1
    case 2
        if strcmp(plot_type,'cross-section')
            disp('Cross-section distributions are only plotted for single products.')
            disp('Plotting first product given.')
            products=products(1);
        end
    otherwise
        switch plot_type
            case 'polar'
                disp('Polar distributions are plotted for at most two products.')
                disp('Plotting two first products given.')
                products=products(1:2);
                n_products=2;
            case 'cross-section'
                disp('Cross-section distributions are only plotted for single products.')
                disp('Plotting first product given.')
                products=products(1);
                n_products=1;
        end
end

% If no axes is given, initialize figure.
if isempty(ax) 
    ax=irf_plot(1);
end

% Plot distributions
switch plot_type
    case 'polar' % Plot polar distribution
        X=cell(1,2);Y=cell(1,2);C=cell(1,2);
        if 0 % For colorbar labels
        product1=products{1};
        product2=products{2};
        if any(strfind(product1,'PSD')); distr1='PSD';
        elseif any(strfind(product1,'DEFlux')); distr1='DEFlux';
        elseif any(strfind(product1,'DPFlux')); distr1='DPFlux';
        end
        if any(strfind(product2,'PSD')); distr2='PSD';
        elseif any(strfind(product2,'DEFlux')); distr2='DEFlux';
        elseif any(strfind(product2,'DPFlux')); distr2='DPFlux';
        end
        end                
        for k=1:n_products % Prepare matrices for surf
            if any(strfind(products{k},'PSD')); distr{k}='PSD';
            elseif any(strfind(products{k},'DEFlux')); distr{k}='DEFlux';
            elseif any(strfind(products{k},'DPFlux')); distr{k}='DPFlux';
            end

            % Preparing different sets of data
            if any(strfind(products{k},'PITCH'))
                [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',products{k}]);

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
                [~,ind]=irf_tlim(tt,tint);

                specrec.f=theta;
                specrec.en=flipdim(en,2); % stigande
                specrec.t=tt;
                specrec.p_label{k}=['Log ' distr ' [' Data_units ']'];
                specrec.f_label='Pitch angle';
                specrec.data=double(data(ind,:,:));
                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN;
                
                r=double(specrec.en)'; % In eV
                rlog=log10(r);
                if isempty(emin)
                    emin= rlog(1);
                end

                %rlog = log10(r-emin); % logarithmic energy levels in eV
                %eminlog=log10(emin);

                r_man = rlog-emin;
                r0=emin;
                theta = double(specrec.f)+90; % turn so that pitch angle zero is on top
                
                % xy-coordinates and c-matrix
                X{k} = r_man*cosd(theta);
                Y{k} = r_man*sind(theta);
                C{k}=flipdim(squeeze(nanmean(log10(specrec.data),1)),2)';
            end
            if any([strfind(products{k},'3DXPH'),...
                    strfind(products{k},'CODIF_HS'),...
                    strfind(products{k},'CODIF_LS'),...
                    strfind(products{k},'HIA_HS_MAG')])
                % Loads data differently depending on product
                if strfind(products{k},'3DXPH'); 
                    res=c_caa_construct_subspin_res_data(['Data__', products{k}]);
                    [caaSEDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', products{k}]);
                    [caaSEDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', products{k}]);
                    SEDL=flipdim(SEDL(1,2:end),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
                    SEDU=flipdim(SEDU(1,2:end),2); en_nan=isnan(SEDU);SEDU(en_nan)=[];                    
                else % Ion products
                    res=c_caa_construct_subspin_res_data(['x3d_ions__', products{k}]);
                    [caaSEDL,~,SEDL]=c_caa_var_get(['delta_plus_energy_table__', products{k}]);
                    [caaSEDU,~,SEDU]=c_caa_var_get(['delta_minus_energy_table__', products{k}]);                            
                    SEDL=flipdim(SEDL(1:end,1),1); en_nan=isnan(SEDL);SEDL(en_nan)=[];
                    SEDU=flipdim(SEDU(1:end,1),1); en_nan=isnan(SEDU);SEDU(en_nan)=[]; 
                end
            
                [~,ind]=irf_tlim(res.tt,tint);
                specrec.t=res.tt;
                specrec.dt=res.dtsampling/2;
                dtheta=(res.theta(2)-res.theta(1))/2;
                specrec.f=[res.theta-dtheta res.theta(end)+dtheta];
                specrec.f_label='Pitch angle';
                specrec.p=res.pitch_angle(ind,:);
                specrec.en=res.en(:);
                specrec.en=[res.en-SEDL' res.en(end)+SEDU(end)'];
                specrec.data=res.data(ind,:,:);
                specrec.p_label{k}=['Log ' distr ' [' res.dataunits ']'];

                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN;

                % Prepare r, and set r0.
                r=double(specrec.en)';
                rlog = log10(r); % energy levels in eV
                if isempty(emin)
                    emin= rlog(1);
                end
                r_man = rlog-emin;
                r0=emin-r_man(1);
                r_man=rlog;

                % Pitch angles, turn so that pitch angle 0 is on top
                theta = double(specrec.f)+90; 

                % xy-coordinates and c-matrix
                X{k} = r_man*cosd(theta); 
                Y{k} = r_man*sind(theta);
                C{k} = squeeze(nanmean(log10(specrec.data),1))';                                           
            end               
        end
        if n_products==1
            X{2}=X{1};Y{2}=Y{1};C{2}=C{1};
        end         
        % Plot data
        surf(ax,X{1},Y{1},X{1}*0,C{1}); hold(ax,'on');
        surf(ax,-flipdim(X{2},2),Y{2},X{2}*0,C{2}); 
        
        view(2); axis(ax,'equal','tight'); shading(ax,'flat'); grid(ax,'off');
        cb=colorbar;
        xlabel('Energy  [eV]')
        ylabel('Energy  [eV]')
        
        t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
        t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
        
        if n_products==1 % Title and labels
            str_product=products{1};
            str_product(strfind(str_product,'_'))=' ';
            title_str={[t1str,'-',t2str,' UT'],str_product};
            ylabel(cb,specrec.p_label{1})
        elseif n_products==2
            str_product1=products{1};
            str_product1(strfind(str_product1,'_'))=' ';
            str_product2=products{2};
            str_product2(strfind(str_product2,'_'))=' ';
            title_str={[t1str,'-',t2str,' UT'],...
                ['Left: ',str_product1],['Right: ', str_product2]};
            ylabel(cb,[specrec.p_label{1} , ' / ', specrec.p_label{2}])
        end
        title(ax,title_str)      
        
        % Ticks
        xticks=log10([1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6])+r0;
        xticks=[xticks(find(xticks>0))];
        xticks=[xticks(find(xticks<r_man(end)))];
        xticklabels=cell(size(xticks));
        for k=1:length(xticklabels)
        xticklabels{k}=num2str(10.^(xticks(k)-r0),'%.1f');
        end         
        xticks=[-flipdim(xticks,2) xticks];
        xticklabels=[flipdim(xticklabels,2) xticklabels];
        yticks=xticks;
        yticklabels=xticklabels; 
        set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
            'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)  
       
        if 1 % Add pitch angle labels
            text(0,rlog(end)-r0-0.5,0,'0^o')
            text(0,-rlog(end)+r0+0.5,0,'180^o')
            text(-rlog(end)+r0+0.5,0,0,'90^o')
            text(rlog(end)-r0-0.5,0,0,'90^o')
        end        
    case 'cross-section' % Plot cross-section distribution
         % Only one product that will be plotted on both sides
           
        product=products{1};

        % Used for colorbar label
        if any(strfind(product,'PSD')); distr='PSD';
        elseif any(strfind(product,'DEFlux')); distr='DEFlux';
        elseif any(strfind(product,'DPFlux')); distr='DPFlux';
        end

        if any(strfind(product,'PITCH'))
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
            phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];

            % Construct subspin reoslution data
            dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            tt=subspintime(dataobject,phi); % Define subspin time from angle phi.
            [~,ind]=irf_tlim(tt,tint);
            % Same for background level
            bgraw=bg.data; bgraw(:,:,:,nan_en)=[];bgraw(:,nan_phi,:,:)=[];
            bgraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            bg=reshape(bgraw,size(bgraw,1)*size(bgraw,2),size(bgraw,3),size(bgraw,4));


            specrec.f=theta;
            specrec.en=flipdim(en,2);
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' Data_units ']'];
            specrec.f_label='Pitch angle';
            specrec.bg=squeeze(nanmean(double(bg(ind,:,:)),1)); % average over time
            specrec.data=squeeze(nanmean(double(data(ind,:,:)),1)); % average over time


            PA0=specrec.data(1,:);
            PA90=nanmean(specrec.data([length(theta)/2 length(theta)/2+1],:),1);
            PA180=specrec.data(end,:);
            PAbg=nanmean(specrec.bg([length(theta)/2 length(theta)/2+1],:),1); % One common zero-count level for all levels

            loglog(ax,specrec.en,PA0,specrec.en,PA90,specrec.en,PA180,specrec.en,PAbg,'--');
            grid(ax,'off');
            %ylabel(cb,specrec.p_label)
            xlabel(ax,'Energy  [eV]')
            ylabel(ax,['Log ', distr ,' ' Data_units ])
            %set(ax,'xlim',[specrec.en(1) specrec.en(end)])

            t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
            prod_str{1}=[t1str,'-',t2str,'UT'];
            prod_str{2}=product;
            text(40,10^7.5,prod_str)   

            legend(ax,'0^o','90^o','180^o','Zero count')      
        end
        if any(strfind(product,'3DXPH'))
            % Load pitch angle data
            res=c_caa_construct_subspin_res_data(['Data__', product]);
            bg=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);

            [~,ind]=irf_tlim(res.tt,tint);
            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            specrec.en=res.en(:);
            specrec.data=res.data(ind,:,:);
            specrec.bg=bg.data(ind,:,:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   

            specrec.data=squeeze(nanmean(specrec.data,1));
            specrec.bg=squeeze(nanmean(specrec.bg,1));
            PA0=specrec.data(1,:)';
            PA90=nanmean(specrec.data([length(specrec.f)/2 length(specrec.f)/2+1],:),1)';
            PA180=specrec.data(end,:)';
            PAbg=nanmean(specrec.bg([length(specrec.f)/2 length(specrec.f)/2+1],:),1)'; % One common zero-count level for all levels

            loglog(ax,specrec.en,PA0,specrec.en,PA90,specrec.en,PA180,specrec.en,PAbg,'--');
            grid(ax,'off');
            %ylabel(cb,specrec.p_label)
            xlabel(ax,'Energy  [eV]')
            ylabel(ax,['Log ', distr ,' ' specrec.p_label ])
            %set(ax,'xlim',[specrec.en(1) specrec.en(end)])

            t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
            prod_str{1}=[t1str,'-',t2str,'UT'];
            prod_str{2}=product;
            text(40,10^7.5,prod_str)   

            legend(ax,'0^o','90^o','180^o','Zero count')  
            energy=specrec.en;
        end
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