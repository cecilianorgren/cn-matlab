function make_pitch(variable_name)
% make pitch angle data from cis/peace polar/azimuthal angles.
% look to c_caa_construct_subspin_res_data.m

%% CODIF_HS does not have pitch angle matrix data, therefore rebinning
% necessary
% CIS sector angels are given in ISR2 reference frame and show particle arrival direction!
theta=11.25:22.5:180; % pitch angles to which rebin
[variable,dataobject,varmat,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
if isempty(variable)
    error('Cannot load the requested CAA variable')
end
phivar=c_caa_var_get(variable.DEPEND_2);phi=phivar.data(1,:);
polarvar=c_caa_var_get(variable.DEPEND_1);polar=polarvar.data(1,:);
envar=c_caa_var_get(variable.DEPEND_3);en=envar.data(1,:);
enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
enlabel=enunits;
% variable.data indices are in order 1) time 2) polar 3) azimuthal 4) energy
pitchangle=variable.data; 
t=varmat.data(:,1);
[tt]=subspintime(dataobject,phi);

dd=regexp(variable_name, '__', 'split');
ic=str2double(dd{2}(2));
c_eval('caa_load C?_CP_FGM_FULL;',ic);
c_eval('B=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');',ic);
B_ISR2=irf_resamp(c_coord_trans('GSE','ISR2',B,'cl_id',ic),tt,'nearest'); % sample to t
b_ISR2=irf_norm(B_ISR2);
cosphi=cos(phi/180*pi);sinphi=sin(phi/180*pi);
cospolar=cos(polar/180*pi);sinpolar=sin(polar/180*pi);
for jphi=1:length(phi),
    for jpolar=1:length(polar),
        % nn vector of given sector in SR2 ref rframe
        nparticle=[cospolar(jpolar).*cosphi(jphi) cospolar(jpolar).*sinphi(jphi) sinpolar(jpolar)];
        pitchsector=acos(dot(b_ISR2(jphi:length(phi):end,2:4),repmat(nparticle,[length(t) 1]),2))*180/pi;
        pitchangle(:,jpolar,jphi,:)=repmat(pitchsector,[1 1 1 length(en)]);
    end
end
variable=fillval_to_nan(variable); % FILLVALs put to NaN
dataraw=ftheta(variable.data,pitchangle,theta);
dataraw=permute(dataraw,[1 2 4 3]); % permute in order [time, pitch, energery, azimuth]

pitch_spectra = squeeze(nansum(dataraw,4));
nEnergies = size(pitch_spectra,3);

specrec.p = pitch_spectra;
specrec.f = theta;
specrec.t = t;

nPanels = 20;
h = irf_plot(nPanels,'newfigure');
%irf_spectrogram(specrec);
if 1
    %%
for oo = 1:nPanels
    hca = h(oo);        
    pcolor(hca,t,theta,log(squeeze(pitch_spectra(:,:,nEnergies-oo)')));    
    hc(oo) = colorbar('peer',hca);
    shading(hca,'flat')
    %caxis(hca,[4.1 6.4]);
    irf_colormap('space')
    irf_legend(hca,[num2str()]
end
end
end

%----------------------------------------------------------
function ftheta=ftheta(fpol,thetapol,theta)
% rebinning to theta vector
% assuming fpol second dimension is polar angle
% rebinnning fpol to ftheta
% assuming thetapol exactly the same size as fpol
ftheta_dim=size(fpol);
ftheta_dim(2)=length(theta);
ftheta=zeros(ftheta_dim);
thetahalfstep=(theta(2)-theta(1))/2;
thetamin=theta-thetahalfstep;
thetamax=theta+thetahalfstep;
for j=1:length(theta),
    ind=(thetapol>thetamin(j)) & (thetapol<thetamax(j));
    fpoltemp=fpol.*ind;
    ind(isnan(fpoltemp))=0;
    fpoltemp(isnan(fpoltemp))=0;
    switch length(ftheta_dim)
        case 3
            ftheta(:,j,:)=sum(fpoltemp,2)./sum(ind,2);
        case 4
            ftheta(:,j,:,:)=sum(fpoltemp,2)./sum(ind,2);
    end
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