function time_str = gcf_time()
ud = get(gcf,'userdata');
%interval_length = diff(ud.subplot_handles(1).XLim);
%interval_start_epochtt = irf_time(ud.t_start_epoch+ud.subplot_handles(1).XLim(1),'epoch>epochtt'); % epochtt
if isfield(ud,'tlim_mva') % irf_minvar_gui
  interval_start_utc = irf_time(ud.tlim_mva(1),'epoch>utc_yyyymmdd_HHMMSSmmm');
  interval_stop_utc = irf_time(ud.tlim_mva(2),'epoch>utc_HHMMSSmmm');
elseif isfield(ud,'ref_satellite') % irf_4_v_gui (_loc)
  interval_start_utc = irf_time(ud.t_start_epoch + ud.h(1).XLim(1),'epoch>utc_yyyymmdd_HHMMSSmmm');
  interval_stop_utc = irf_time(ud.t_start_epoch + ud.h(1).XLim(2),'epoch>utc_HHMMSSmmm');
else
  interval_start_utc = irf_time(ud.t_start_epoch + ud.subplot_handles(1).XLim(1),'epoch>utc_yyyymmdd_HHMMSS');
  interval_stop_utc = irf_time(ud.t_start_epoch + +ud.subplot_handles(1).XLim(2),'epoch>utc_HHMMSS');
end
time_str = sprintf('%s_%s',interval_start_utc,interval_stop_utc);

