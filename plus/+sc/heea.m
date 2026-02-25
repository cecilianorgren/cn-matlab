function angle = heea(B,sc)
% SC.HEEA Calculate angle between peace detectors and magnetic field.

% take out time interval for which we need to download phase data.
tint = torow(B([1 end],1));

% download phase data
c_eval('[tt,phase_data] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], tint(1), diff(tint));',sc);

% resample it so it has a "normal" time difference
phase=c_phase(B(:,1),[tt phase_data]);

% resample it to input
phase = irf_resamp(phase,B); % phase time series resampled to B

% get heea and leea phases
phase_heea=data.phase/180*pi-(30)/180*pi;
phase_leea=phase_heea+pi;

% get detector (normalized) vectors in ds reference frame
rheea=irf_norm([cos(phase_heea) sin(phase_heea);cos(phase_heea-dphi) sin(phase_heea-dphi);cos(phase_heea+dphi) sin(phase_heea+dphi)]);
rleea=irf_norm([cos(phase_leea) sin(phase_leea);cos(phase_leea-dphi) sin(phase_leea-dphi);cos(phase_leea+dphi) sin(phase_leea+dphi)]);

% transform to isr2
c_coord_trans('dsc','isr2',rheea,'cl_id',sc);
c_coord_trans('dsc','isr2',rleea,'cl_id',sc);
        
