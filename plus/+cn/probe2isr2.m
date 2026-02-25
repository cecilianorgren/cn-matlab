function out = probe2isr2(in,sc,probe)
% CN.PROBE2ISR2     Transforms measurments made by probes in spinning frame
%                   to ISR2 frame. Used to derive electric field from raw 
%                   data.
%
%   output = CN.PROBE2ISR2(input,spacecraft,probe);
%       input - product thats in the spinning spacecraft system
%       output - product, transformed to ISR2 system
%       sc - spacecraft
%       probe - which probe the input is measured along


% take out time interval for which we need to download phase data.
tint = torow(in([1 end],1));

% download phase data
c_eval('[tt,phase_data] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], tint(1), diff(tint));',sc);

% resample it so it has a "normal" time difference
phase=c_phase(in(:,1),[tt phase_data]);

% resample it to input
phase = irf_resamp(phase,in);

% calculate phase of individual probes
phase_p1 = irf_tappl(phase,'/180*pi + 3*pi/4');
phase_p3 = irf_tappl(phase_p1,' - pi/2');
phase_p2 = irf_tappl(phase_p1,' + pi');
phase_p4 = irf_tappl(phase_p1,' + pi/2');


% time series of in input in dsc coordinate system, which is same as sr2
c_eval('isr_out = [phase_p?(:,1) in(:,2).*cos(phase_p?(:,2)) in(:,2).*sin(phase_p?(:,2)) phase_p?(:,2)*0];',probe);

% transforming input to isr2 system, not dsc=sr2 and dsi=isr2
c_eval('out=c_coord_trans(''dsc'',''isr2'',isr_out,''cl_id'',sc);',probe);
