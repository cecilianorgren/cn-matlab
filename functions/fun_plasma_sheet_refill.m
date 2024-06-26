function refill_time = fun_plasma_sheet_refill(vlb,nlb,nps,Lps)
% refill_time = fun_plasma_sheet_refill(vlb,nlb,nps,zps)
% refill_time in seconds
% refill_time = Lps*nps/(vlb*nlb)
 
flux_lobe = (vlb*1e3)*(nlb*1e6); % m/s*m-3
content_ps = (Lps*1e3)*(nps*1e6); % m*m-3

replacement_rate = flux_lobe/content_ps; % (m/s*m-3)/(m*m-3)

refill_time = 1/replacement_rate; % s
%refill_time = refill_time/60/60; % h