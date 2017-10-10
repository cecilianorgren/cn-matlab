
%% Wave analysis
if 0
  %% Field aligned wavelets
  c_eval('facE? = irf_convert_fac(gseE?,gseB?,[1 0 0]);',ic)
  c_eval('facB?scm = irf_convert_fac(gseB?scm,gseB?,[1 0 0]);',ic)
  
  c_eval('wavE?fac = irf_wavelet(facE?.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavE?fac.f_units = ''Hz''; wavE?fac.f_label = ''f [Hz]''; wavE?fac.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
  c_eval('wavB?fac = irf_wavelet(facB?scm.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavB?fac.f_units = ''nT''; wavB?fac.f_label = ''f [Hz]''; wavB?fac.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
%%
if 0
  ic = 1:4;
  tic
  c_eval('wavE? = irf_wavelet(gseE?.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
  c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
  c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
    toc
%%
tintPol = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:52.00Z');
c_eval('polarization? = irf_ebsp(gseE?.tlim(tintPol),gseB?scm.tlim(tintPol),gseB?.tlim(tintPol),gseB?.tlim(tintPol),gseR?brsttime.tlim(tintPol),[10 2200],''polarization'',''fac'');',1)
end
end