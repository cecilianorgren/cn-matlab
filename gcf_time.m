function time_str = gcf_time()
ud = get(gcf,'userdata');
%interval_length = diff(ud.subplot_handles(1).XLim);
%interval_start_epochtt = irf_time(ud.t_start_epoch+ud.subplot_handles(1).XLim(1),'epoch>epochtt'); % epochtt
interval_start_utc = irf_time(ud.t_start_epoch + ud.subplot_handles(1).XLim(1),'epoch>utc_yyyymmdd_HHMMSS');
interval_stop_utc = irf_time(ud.t_start_epoch + +ud.subplot_handles(1).XLim(2),'epoch>utc_HHMMSS');
time_str = sprintf('%s_%s',interval_start_utc,interval_stop_utc);

