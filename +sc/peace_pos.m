function out = peace_pos(B,sc,varargin)
% SC.PEACE_POS Calculate angle between peace detectors and magnetic field.
%   SC.PEACE_POS(B,sc,detector)
%          B - magnetic field in ISR2 coordinates
%         sc - spacecraft number, needed for coordinate tranfsormation
%   detector - specify output, can be 'heea' [time angle], 
%                                     'leea' [time angle]
%                                     'both' [time angle_heea angle_leea]

% take out time interval for which we need to download phase data.
tint = torow(B([1 end],1));

% download phase data
c_eval('[tt,phase_data] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], tint(1), diff(tint));',sc);

% resample it so it has a "normal" time difference
phase=c_phase(B(:,1),[tt phase_data]);

% resample it to input
phase = irf_resamp(phase,B); % phase time series resampled to B

% get heea and leea phases
phase_heea=phase(:,2)/180*pi-(30)/180*pi;
phase_leea=phase_heea+pi;

% get detector (normalized) vectors in ds reference frame
rheea=[irf_norm([cos(phase_heea) sin(phase_heea)])]; % only in plance component
rleea=[irf_norm([cos(phase_leea) sin(phase_leea)])]; % only in plance component

% transform B to isr and normalize, its easier that way, only need angle anyway
dB = irf_norm(c_coord_trans('isr2','dsc',B,'cl_id',sc)); 
dB_plane_norm = sqrt(dB(:,2).^2+dB(:,3).^2);
       
% get angles between, heea/leea and B, in spin plane, and out of plane 
plane_angle_leea = acosd(dot(dB(:,[2 3]),rleea,2)./dB_plane_norm);
plane_angle_heea = acosd(dot(dB(:,[2 3]),rheea,2)./dB_plane_norm);

% prepare output
if nargin>2
    detector=varargin{1};
    switch lower(detector)
        case 'heea'            
            out = [B(:,1) plane_angle_heea];
        case 'leea'            
            out = [B(:,1) plane_angle_heea];
        case 'both'
            out = [B(:,1) plane_angle_heea plane_angle_leea];
    end
end 