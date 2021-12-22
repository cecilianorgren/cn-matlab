function out  = edi_event_manual_dt_phi
% New wave properties
manual = edi_event_manual_dt;

t_refs = EpochTT(cat(1,manual.t_ref));

tint = t_refs([1 end]) + [-10 10];
t0 = t_refs(1);

ic = 1:4;
doLoad = 1;
if doLoad
  c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',ic);
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)  
end

c_eval('gseR? = gseR?.resample(t0);',ic);
R0 = (gseR1.data + gseR2.data + gseR3.data + gseR4.data)/4; % center of tetrahedron
c_eval('gseR?rel = gseR? - R0;',ic); % distance in km from center of tetrahedron

%% Get remaining ESW properties: length scale, potential
irf_log OFF

data = [];
colors = mms_colors('1234');
% Run through to estimate peak-to-peak length scales, and possible
% potentials
for ieh = 1:numel(manual)
  T = 0.006;
  tref = EpochTT(manual(ieh).t_ref);
  dt = manual(ieh).dt;
  tminus = EpochTT(manual(ieh).tminus);
  tplus = EpochTT(manual(ieh).tplus);
  tp1 = tref + manual(ieh).tp1;
  tp2 = tref + manual(ieh).tp2;
  tpp = tp2-tp1;
  lmn = [manual(ieh).perp1; manual(ieh).perp2; manual(ieh).par]; 
  gsev = manual(ieh).v; % velocity is given in gse coordinate system
  pppv = gsev*lmn';
  pppvnorm = pppv/norm(pppv); % unit vector of velocity
  vpar = pppv(3);  
  tint_eh_all = tref + [min(dt) max(dt)]+0.5*T*[-1 1];
  c_eval('E? = gseE?.tlim(tint_eh_all)*lmn'';',1:4) % use same time interval for all
  c_eval('R?rel = gseR?rel.data*lmn'';',1:4)
  gseBav = manual(ieh).gseBav;
  pppB = gseBav*lmn';
  
  c_eval('intE? = irf_integrate(E?.z,tref+dt(?));',1:4) % (mV/m)*s
  c_eval('phi? = vpar*intE?;',1:4) % (km/s)*(mV/m)*s = V  
  
  % Try to estimate what the potential should be based on the value at
  % tminus and tplus  
  %c_eval('tminus? = tminus + dt(?);',1:4)
  %c_eval('tplus?  = tplus  + dt(?);',1:4)    
  c_eval('phiminus? = phi?.resample(tminus + dt(?));',1:4)
  c_eval('phiplus?  = phi?.resample(tplus + dt(?));',1:4)
  c_eval('phiminusplus?  = irf.ts_scalar([tminus tplus] + dt(?),[phiminus?.data,phiplus?.data]);',1:4)  
  % c_eval('phiminusplus?  = phi?.resample([tminus tplus] + dt(?));',1:4)
  % This above is not working, and I have no idea why...
  c_eval('phimean(?) = 0.5*(phiplus?.data + phiminus?.data);',1:4)  
  
  data(ieh).tpp = tpp;  
  data(ieh).lpp = tpp*vpar;
  data(ieh).phiminus = [phiminus1.data phiminus2.data phiminus3.data phiminus4.data];
  data(ieh).phiplus = [phiplus1.data phiplus2.data phiplus3.data phiplus4.data];  
  data(ieh).phimean = phimean;
  data(ieh).vpar = vpar;
  data(ieh).pitchangle = acosd(pppvnorm(3));
  
  
  % Plot
  if 0
    figure(79)
    h = irf_plot(3);

    hca = irf_panel('Epar');  
    set(hca,'ColorOrder',mms_colors('1234'))
    comp = 'z';
    irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp')
    hca.YLabel.String = {'E_{||} (mV/m)'};

    hca = irf_panel('Epar dt');  
    set(hca,'ColorOrder',mms_colors('1234'))
    comp = 'z';
    irf_plot(hca,{E1.(comp),E2.(comp),E3.(comp),E4.(comp)},'comp','dt',dt)
    hca.YLabel.String = {'E_{||} (mV/m)'};
    irf_legend(hca,{sprintf('<tpp>=%.0f',tpp(1));...
                    sprintf('<tpp>=%.0f',tpp(2));...
                    sprintf('<tpp>=%.0f',tpp(3));...
                    sprintf('<tpp>=%.0f',tpp(4))},[0.02 0.02])

    hca = irf_panel('phi');    
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt)
    %set(hca,'ColorOrder',mms_colors('11223344'))
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{phiminusplus1,phiminusplus2,phiminusplus3,phiminusplus4},'comp','dt',dt);
    irf_legend(hca,{sprintf('<phi1>=%.0f',phimean(1));...
                    sprintf('<phi2>=%.0f',phimean(2));...
                    sprintf('<phi3>=%.0f',phimean(3));...
                    sprintf('<phi4>=%.0f',phimean(4))},[0.02 0.02])
     %c_eval('hlines(?).Marker = ''o'';',1:2)
  %   set(hca,'ColorOrder',mms_colors('22'))
  %   hlines = irf_plot(hca,{phiminus2,phiplus2},'comp');
  % %   c_eval('hlines(?).Marker = ''o'';',1:2)   
  %   hlines = irf_plot(hca,{phiminus3,phiplus3},'comp');
  %  %  c_eval('hline.s(?).Marker = ''o'';',1:2)   
  %   hlines = irf_plot(hca,{phiminus4,phiplus4},'comp');
  %   % c_eval('hlines(?).Marker = ''o'';',1:2)
  %   hca.YLabel.String = {'\phi (V)'};

    h(1).Title.String = sprintf('ieh = %g',ieh);
    irf_zoom(h,'x',tint_eh_all)
    irf_pl_mark(h,tref,[0 0 0])
    irf_pl_mark(h,[tminus tplus],[1 1 0])
    c_eval('irf_pl_mark(h,tp1(?),colors(?,:))',1:4)
    c_eval('irf_pl_mark(h,tp2(?),colors(?,:))',1:4)

    fprintf('ieh = %g\n',ieh)

    %[tphi,levphi] = get_time(2,'epochtt');
    %fprintf('manual(ieh).tminus = ''%s'', manual(ieh).tplus = ''%s'' \n',tphi(1).utc,tphi(2).utc)

    %[tpeaks,levpeaks] = get_time(8,'epochtt');
    %tp1 = [tpeaks(1),tpeaks(3),tpeaks(5),tpeaks(7)]-tref;
    %tp2 = [tpeaks(2),tpeaks(4),tpeaks(6),tpeaks(8)]-tref;  
    %fprintf('manual(ieh).tp1 = [%8.5f,%8.5f,%8.5f,%8.5f];\n',tp1(1),tp1(2),tp1(3),tp1(4))
    %fprintf('manual(ieh).tp2 = [%8.5f,%8.5f,%8.5f,%8.5f];\n',tp2(1),tp2(2),tp2(3),tp2(4))
    %pause
  end
end
irf_log ON
out = data;