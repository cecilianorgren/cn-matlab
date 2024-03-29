%% Outside loop
% Everything that only needs to be done once
localuser = datastore('local','user');
units = irf_units;

% Load MMS data
if 0 % load from .mat file
  load /Users/cecilia/Data/20170706_135303_basic_eh
elseif 1 % load data
  ic = 1:4;
  tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
  mms.db_init('loca_file_db','/Volumes/Nexus/data');x
  db_info = datastore('mms_db');   
  c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
  c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
  c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
  c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
  c_eval('intEdt? = irf_integrate(gseE?par);');    
  mms.load_data_edi;  
  c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
  % Make reduced distribution
  tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.617Z');
  tintZoom = tint_phi + [-2 2];   
  strTintZoom = [irf_time(tintZoom(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tintZoom(2),'epochtt>utc_HHMMSS')];
  eint = [000 40000];
  vint = [-Inf Inf];
  eDist = ePDist1.tlim(tintZoom).elim(eint);   
  scpot = scPot1.resample(eDist);
  scpot_margin = 1.0; % keep in mind that this also affects the velocity at lower energies
  lowerelim = scpot/scpot*50;
  eLine = dmpaB1.resample(eDist).norm;
  tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B  
end

tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.620Z');
t0 = tint_phi(1);

% FPI distribution
ts_f_fpi = ef1D.tlim(tint_phi);
v_fpi = mean(ts_f_fpi.depend{1},1);
f_fpi = ts_f_fpi.data;

% EDI flux
edi_nodes = 1;
c_eval('ts_edi_flux180_mms? = irf.ts_scalar(flux180_mms?.tlim(tint_phi).time,mean(flux180_mms?.tlim(tint_phi).data(:,edi_nodes),2));')
c_eval('ts_edi_flux0_mms? = irf.ts_scalar(flux0_mms?.tlim(tint_phi).time,mean(flux0_mms?.tlim(tint_phi).data(:,edi_nodes),2));')

% EDI energy and corresponding velocity
E_edi = 500; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s

% Prepare electric field
% later we need just to multiply with vph
ffilt = 30; % Hz
c_eval('Etoint? = gseE?par;')
c_eval('Etoint? = Etoint?.filt(ffilt,0,[],3).tlim(tint_phi);');    
c_eval('intEdt? = irf_integrate(Etoint?);');    
minpeakdistance = 30;
c_eval('[PKS?,LOCS?,W?] = findpeaks(intEdt?.data,''MinPeakDistance'',minpeakdistance);')
c_eval('intEdt?_detrend = intEdt?; intEdt?_detrend.data = detrend(intEdt?_detrend.data,''linear'',LOCS?);')

% Plotting options
doT = 1; % otherwise plot x;
x_vec = intEdt1.time - t0; % seconds
dx_vec = x_vec(2)-x_vec(1);
nx = numel(x_vec);
x_vec_diff1 = x_vec(1:end-1)+0.5*dx_vec;
x_vec_diff2 = x_vec(2:end-1);
        
% Remove phi baselevel
%c_eval('phi_baselevel = interp_linear_piecewise(intEdt?.data,x_vec,x_vec(LOCS?));')
c_eval('ts_locs? = irf.ts_scalar(intEdt?.time([1; LOCS?; end]),intEdt?.data([1; LOCS?; end]));')
%c_eval('phi_baselevel? = interp_linear_piecewise(intEdt?.data([1; LOCS?; end]),x_vec([1; LOCS?; end]),x_vec);')
%c_eval('ts_phi_baselevel? = irf.ts_scalar(intEdt?.time,phi_baselevel?);')
c_eval('ts_phi_baselevel? = ts_locs?.resample(intEdt?);')
c_eval('intEdt?_detrend = intEdt?-ts_phi_baselevel?;')

% EH properties
% observed/measured properties (konrad)
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
% charge separation assuming potential structure is single gaussian
dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc


% velocity grid
vmax = 50000e3;
nv = 3000;
v_vec = linspace(-vmax,vmax,nv);
dv = v_vec(2) - v_vec(1);

% flux at edi energy/velocity interval
vind_edi_0 = intersect(find(v_vec>v_edi_minus),find(v_vec<v_edi_plus));
vind_edi_180 = intersect(find(v_vec<-v_edi_minus),find(v_vec>-v_edi_plus));

% How to use phase velocity
option_vph = 'interp linear piecewise';
option_vph = 'constant';

%% Inside loop
doPrint = 0;

% Phase speed and potential zero level shift
vph = -9000e3; % m/s, representative phase velocity
phi_shift = 00; % to keep potential > 0
phi_scaling = 1.0; % if we underestimate phi

% Parameters to vary
% two populations, f1 and f2
nPop = 2;
ntot = 0.035*1e6;
R_min = 0.50; R_max = 0.7; nR = 3; vec_R = linspace(R_min,R_max,nR); % rato of f1 to f2
vec_N1 = ntot*vec_R; 
vec_N2 = ntot*(1-vec_R);
T1_min = 200; T1_max = 300; nT1 = 3; vec_T1 = linspace(T1_min,T1_max,nT1);
%T2_min = 50; T2_max = 400; nT2 = 10; vec_T2 = linspace(T2_min,T2_max,nT2);
T2_min = 4000; T2_max = 2000; nT2 = 1; vec_T2 = linspace(T2_min,T2_max,nT2);
Vd1_min = -10000*1e3; Vd1_max = -7000*1e3; nVd1 = 3; vec_Vd1 = linspace(Vd1_min,Vd1_max,nVd1);
%Vd2_min = 0000*1e3; Vd2_max = 4000*1e3; nVd2 = 3; vec_Vd2 = linspace(Vd2_min,Vd2_max,nVd);
Vd2_min = 4000*1e3; Vd2_max = 3000*1e3; nVd2 = 1; vec_Vd2 = linspace(Vd2_min,Vd2_max,nVd2);
Beta_min = -0.80; Beta_max = -0.2; nBeta = 5; vec_Beta = linspace(Beta_min,Beta_max,nBeta);

mms_id = 1; % chose spacecraft to compare to
c_eval('LOCS = LOCS?; PKS = PKS?;',mms_id)
c_eval('obs_eh_xvec = obs_t0_epoch_mms?-t0;',mms_id)
c_eval('ts_edi_flux0 = ts_edi_flux0_mms?;',mms_id)
c_eval('ts_edi_flux180 = ts_edi_flux180_mms?;',mms_id)
nN = nR*nT1*nT2*nVd1*nVd2*nBeta;

for iR = 1%:nR
  for iT1 = 1%:nT1
    for iT2 = 1%:nT2
      for iVd1 = 1%:nVd1
        for iVd2 = 1%:nVd2
          for iBeta = 1%:nBeta
            % Plasma properties
            R = vec_R(iR);
            n = ntot*[R 1-R];
            T = [vec_T1(iT1) vec_T2(iT2)]; % eV
            vd = [vec_Vd1(iVd1) vec_Vd2(iVd2)]; % m/s
            beta = vec_Beta(iBeta); % in principle, this could be varied between the two species
            

            % Physical parameters
            vt = sqrt(2*units.e*T./units.me); % m/s

            % Potential from observed E
            c_eval('phi_timeline = intEdt?_detrend.time;',mms_id)
            c_eval('phi?_detrend = intEdt?_detrend*vph*1e-3*phi_scaling;',mms_id)
            c_eval('phi?_detrend_shift = phi?_detrend + phi_shift;',mms_id)
            c_eval('phi?_detrend_shift.data(phi?_detrend_shift.data<0) = 0;',mms_id)                
            c_eval('phi_vec = phi?_detrend_shift.data;',mms_id)

            c_eval('epar_vec = Etoint?.data;',mms_id)
            c_eval('epar_vec_nofilt = gseE?par.tlim(tint_phi).data;',mms_id)

            switch option_vph
              case 'gradual'
                vph_vec = vph*(1+0.5*x_vec/max(x_vec));
              case 'constant'
                vph_vec = repmat(vph,1,nx);
              case 'interp linear piecewise' % piecewise linear interpolation of vph
                c_eval('tmp_time = obs_t0_epoch_mms?-t0;',mms_id)
                tmp_time = torow([tint_phi(1)-t0; tmp_time]);
                c_eval('tmp_data = obs_vph?.data*1e3;',mms_id) % m/s
                tmp_data = torow([vph; tmp_data]);        
                vph_vec = interp_linear_piecewise(tmp_data,tmp_time,x_vec);
                vph_vec = smooth(vph_vec,numel(x_vec)/50);
                %plot(x_vec,vph_vec,'.',tmp_time,tmp_data,'*')
            end

            % adjust phi incase vph is variable (vph is used to get phi: phi = eint*vph)
            phi_vec = phi_vec.*reshape(vph_vec,size(phi_vec))/vph;

            t_to_x = -vph;
            if doT  
              dx_vec = x_vec(2) - x_vec(1);
              dx = dx_vec*t_to_x;
            else
              dx_vec = x_vec(2) - x_vec(1);
              dx = x_vec(2) - x_vec(1);
            end

            [X,V] = meshgrid(x_vec,v_vec);
            VPH = repmat(torow(vph_vec),nv,1);
            PHI = repmat(torow(phi_vec),nv,1);
                       
            % charge density from observed phi        
            obs_density_diff = -diff(epar_vec*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
            obs_density_diff_nofilt = -diff(epar_vec_nofilt*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
            dn = obs_density_diff; dn = [0; dn];
            
            if 0 % Abel, Bernstein1957
              [Fflat,Fflat_free,Fflat_trap] = get_f_flat(V',n,vt,vd,PHI',VPH');
              nfree = nansum(Fflat_free,2)*dv;
              ntrap_flat = nansum(Fflat_trap,2)*dv;
              ntrap = ntot - torow(nfree) + torow(dn);
              dntrap = ntrap - torow(ntrap_flat);
              [Fabel,Fabel_free,Fabel_trap] = get_f_abel(V',n,vt,vd,PHI',VPH',dntrap);
              Fabel = Fabel + Fflat_trap;
              F = Fabel';
              ifree = find(Fflat==0);
              itrap = find(Fflat>0);
            else % Schamel
              % model phase space density distribution                      
              % F = fe_schamel(V,n,vt,vd,PHI,VPH,beta); % distribution function, from Schamel 1986                        
              F = zeros(size(V));
              Ffree = zeros(size(V));
              Ftrap = zeros(size(V));
              Fmax = zeros(size(V));
              for iPop = 1:nPop
                Fmax_tmp = fe_maxwellian(V,n(iPop),vt(iPop),vd(iPop)); % distribution function, from Schamel 1986
                c_eval('Fmax? = Fmax_tmp;',iPop)
                Fmax = Fmax + Fmax_tmp;

                [Ftmp,Ftmp_all] = fe_schamel(V,n(iPop),vt(iPop),vd(iPop),PHI,VPH,beta); % distribution function, from Schamel 1986
                c_eval('F? = Ftmp;',iPop)              
                F = F + Ftmp; 
                Ffree = Ffree + Ftmp_all.free; 
                Ftrap = Ftrap + Ftmp_all.trap; 
              end                    
            end
            
            %% 'moments'
            Fdv = F*dv;
            FV = F.*V; % (s1/m4)*(m1/s1) = (1/m3)  
            FVdv = FV*dv; % (s1/m4)*(m1/s1)^2 = (1/m2s)
            %Fdv_trap = F*dv; Fdv_trap(Ftmp_all.ifree) = NaN; % (s1/m4)*(m1/s1)^2 = (1/m2s)
            %Fdv_free = F*dv; Fdv_free(Ftmp_all.itrap) = NaN; % (s1/m4)*(m1/s1)^2 = (1/m2s)
            Fdv_trap = F*dv; Fdv_trap(ifree) = NaN; % (s1/m4)*(m1/s1)^2 = (1/m2s)
            Fdv_free = F*dv; Fdv_free(itrap) = NaN; % (s1/m4)*(m1/s1)^2 = (1/m2s)

            
            % 'integrals'
            sumFdv = nansum(Fdv,1); % n: (s1/m4)*(m1/s1) = (1/m3) 
            sumFVdv = nansum(FVdv,1); % nv: (s1/m4)*(m1/s1)*(m1/s1) = (1/m2s)        
            mod_density = sumFdv;
            mod_density_free = nansum(Fdv_free,1);
            mod_density_trap = nansum(Fdv_trap,1);
            mod_velocity = sumFVdv./sumFdv;
            ts_mod_density = irf.ts_scalar(phi_timeline,mod_density);
      
            % 'averages'
            mod_f_average = mean(F,2);
            mod_fmax_average = mean(Fmax,2);
            mod_density_average = mean(mod_density);
            mod_velocity_average = mean(mod_velocity);

            % charge density from observed phi        
            %obs_density_diff = -diff(epar_vec*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
            %obs_density_diff_nofilt = -diff(epar_vec_nofilt*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
            obs_density = mod_density_average + obs_density_diff; % ni + ne - ni = ne,  assuming ion density is unperturbed  
            ts_obs_density = irf.ts_scalar(phi_timeline(1:end-1) + 0.5*(phi_timeline(2)-phi_timeline(1)),obs_density);
               
            % get net from phi
            nef = mod_density_free;
            n_diff =  torow([obs_density_diff(1); obs_density_diff]);
            net = ntot - nef + n_diff;
            n_diff_nofilt =  torow([obs_density_diff_nofilt(1); obs_density_diff_nofilt]);
            net_nofilt = ntot - nef + n_diff_nofilt;
            %net_di0ff_ephi = diff(net)/units.e

            % BGK approach
            [Ftmp_bgk,Ftmp_bgk_all] = fe_bgk(V,n,vt,vd,PHI,VPH,obs_density_diff,phi_vec); % distribution function, from Schamel 1986
                        
            % model flux
            model_flux_edi_0 = nansum(FVdv(vind_edi_0,:));
            model_flux_edi_180 = nansum(FVdv(vind_edi_180,:));
            units_scale = 1e-4; % 1/m2 -> 1/cm2
            ts_mod_flux_0 = irf.ts_scalar(phi_timeline,model_flux_edi_0*units_scale);
            ts_mod_flux_180 = irf.ts_scalar(phi_timeline,abs(model_flux_edi_180)*units_scale);

            % observed flux
            % ts_edi_flux0 
            % ts_edi_flux180

            % Decide quality of comparison
            % density,
            [Q_Flux_xcorr,Q_Flux_area] = quality_of_fit(ts_edi_flux180.data,ts_mod_flux_180.resample(ts_edi_flux180).data);
            [Q_Density_xcorr,Q_Density_area] = quality_of_fit(ts_obs_density.data,ts_mod_density.resample(ts_obs_density).data);
            
            if 1 % plot, for diagnostics etc
              %% time series
              fig = figure(33);
              scrsz = get(groot,'ScreenSize');
              figurePostition = scrsz; 
              figurePostition(1) = 100;
              figurePostition(2) = 100;
              figurePostition(3) = figurePostition(3)*0.8; 
              figurePostition(4) = figurePostition(4)*0.8;
              fig.Position = figurePostition;
              clear h;
      
              nrows = 9;
              ncols = 2;
              npanels = nrows*1;
              isub = 1;
              for irow = 1:nrows
                for icol = 1
                  ipanel = 1+(irow-1)*2;
                  h(isub) = subplot(nrows,ncols,ipanel);
                  isub = isub + 1;
                end
              end
              
              isub = 1;

              vlim_f = 30;
              if 1 % Epar, plot
                hca = h(isub); isub = isub + 1;
                plot(hca,x_vec,epar_vec);  
                hca.YLabel.String = {'E_{||}','(mV/m)'};
              end
              if 1 % PHI, plot
                hca = h(isub); isub = isub + 1;
                plot(hca,x_vec,phi_vec);  
                if 1
                  hold(hca,'on')
                  plot(hca,x_vec(LOCS),phi_vec(LOCS),'*')
                  hold(hca,'off')
                end
                irf_legend(hca,{sprintf('v_{ph}=%gx10^3 km/s, phi_{shift}=%g V, f_{filt,E}=%g Hz',vph*1e-6,phi_shift,ffilt)},[0.01 0.99],'color',[0 0 0]);
                hca.YLabel.String = {'\phi','(V)'};  
              end
              if 1 % diff E (Poisson) (density)
                hca = h(isub); isub = isub + 1;  
                if 1
                  nscale = 1e-3;
                  plot(hca,x_vec_diff1,obs_density_diff*1e-6/nscale); % 1e-6 from 1/m3 > 1/cm3
                  %legend(hca,{'n_{obs,Poisson,1D}'},'location','eastoutside');
                end
                hca.YLabel.String = {'\delta n',sprintf('(10^{%.0f} cm^{-3})',log10(nscale))};
              end
              if 1 % F
                hca = h(isub); isub = isub + 1;
                pcolor(hca,X,V*1e-6,(F));
                shading(hca,'flat')
                hca.YLabel.String = {'v','(10^3 km/s)'};
                hcb = colorbar('peer',hca);
                hcb.YLabel.String = 'f (s^1/m^4)';
                colormap(hca,cn.cmap('white_blue'))
                %colormap(hca,cn.cmap('blue_white'))
                hca.YLim = vlim_f*[-1 1];
                
                if 1 % EDI energies
                  hold(hca,'on')
                  line_color = [0.5 0.5 0.5]; [0.9290    0.6940    0.1250];
                  hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
                  for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end
                  hlines = plot(hca,...
                    x_vec([1 end]),(v_edi_plus)*1e-6*[1 1],...
                    x_vec([1 end]),(v_edi_minus)*1e-6*[1 1],...
                    x_vec([1 end]),(-v_edi_plus)*1e-6*[1 1],...
                    x_vec([1 end]),(-v_edi_minus)*1e-6*[1 1],...
                    'LineWidth',1.5);
                  for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
                  irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);  
                  hold(hca,'off')
                end
                if 1 % model phase velocity
                  hold(hca,'on')
                  line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
                  hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
                  irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
                  hold(hca,'off')
                end
                if 1 % observed phase velocity
                  hold(hca,'on')
                  hlines = plot(hca,obs_eh_xvec,obs_velocity*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
                  irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
                  hold(hca,'off')
                end  
                str_info = {'unperturbed f:';...
                  ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
                  ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
                  ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
                  sprintf('beta_{Schamel}=%g',beta);...
                  };
                irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
                %irf_legend(hca,{sprintf('T_{bg}=%g eV, n_{bg}=%g cc, v_{d,bg}=%g km/s, beta_{Schamel}=%g',T,n*1e-6,vd*1e-3,beta)},[0.99 0.99],'color',hlines(1).Color);    
              end
              if 0 % e psd vpar
                hca = h(isub); isub = isub + 1;
                specrec = ef1D.tlim(tint_phi).specrec('velocity_1D','10^3 km/s');
                time_fred = ef1D.tlim(tint_phi).time-tint_phi(1);
                imagesc(hca,time_fred,specrec.f(1,:),specrec.p');
                hcb = colorbar('peer',hca);
                colormap(cn.cmap('white_blue'))
                hca.CLim = [0 3]*1e-3; 

                hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
                hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
                irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
                irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
              end
              if 0 % sum F (density,free,trap,all)
                hca = h(isub); isub = isub + 1;  
                if 1
                  plot(hca,x_vec([1 end]),ntot*1e-6*[1 1],x_vec,nef*1e-6,x_vec,net*1e-6,x_vec,(ntot+n_diff)*1e-6); % 1e-6 from 1/m3 > 1/cm3
                  legend(hca,{'n_{mod,\phi=0}','n_{mod,free}','n_{mod,trap}'},'location','eastoutside');
                end
                hca.YLabel.String = {'n','(cm^{-3})'};
              end
              if 1 % sum F (density,free,trap,all)
                hca = h(isub); isub = isub + 1;                  
                plot(hca,x_vec([1 end]),ntot*1e-6*[1 1],x_vec,mod_density*1e-6,x_vec,mod_density_free*1e-6,x_vec,mod_density_trap*1e-6,x_vec,net*1e-6); % 1e-6 from 1/m3 > 1/cm3
                legend(hca,{'n_{mod,\phi=0}','n_{mod}','n_{mod,free}','n_{mod,trap}','n_{mod,\phi=0} - n_{mod,free}+\delta n'},'location','eastoutside');
                hca.YLabel.String = {'n','(cm^{-3})'};
              end
              if 1 % sum F (density)
                hca = h(isub); isub = isub + 1;  
                if 1
                  plot(hca,x_vec,mod_density*1e-6,x_vec_diff1,obs_density*1e-6); % 1e-6 from 1/m3 > 1/cm3
                  legend(hca,{'n_{mod}','<n_{mod}> + \delta n'},'location','eastoutside');
                  irf_legend(hca,{sprintf('Q_X=%.2f, Q_A=%g',Q_Density_xcorr,Q_Density_area)},[0.01 0.1],'color',[0 0 0])
                else
                  plot(hca,x_vec([1 end]),ntot*1e-6*[1 1],x_vec,mod_density*1e-6,x_vec([1 end]),mean(mod_density*1e-6)*[1 1],x_vec_diff1,obs_density*1e-6); % 1e-6 from 1/m3 > 1/cm3
                  legend(hca,{'n_{mod,\phi=0}','n_{mod}','<n_{mod}>','n_{obs,Poisson,1D}'},'location','eastoutside');
                end
                hca.YLabel.String = {'n','(cm^{-3})'};
              end
              if 1 % flux: F*v
                hca = h(isub); isub = isub + 1;
                pcolor(hca,X,V*1e-6,FV);
                shading(hca,'flat')
                hca.YLabel.String = {'v','(10^3 km/s)'};
                hcb = colorbar('peer',hca);
                hcb.YLabel.String = 'f*v (1/m^3)'; 
                hca.CLim = hca.CLim(2)*[-1 1]; 
                hcb.YLim = hca.CLim;
                hca.YLim = vlim_f*[-1 1];
                if 1 % EDI energies
                  hold(hca,'on')
                  hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
                  for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
                  hlines = plot(hca,...
                    x_vec([1 end]),(v_edi_plus)*1e-6*[1 1],...
                    x_vec([1 end]),(v_edi_minus)*1e-6*[1 1],...
                    x_vec([1 end]),(-v_edi_plus)*1e-6*[1 1],...
                    x_vec([1 end]),(-v_edi_minus)*1e-6*[1 1],...
                    'LineWidth',1.5);
                  for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
                  irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
                  hold(hca,'off')
                end
                if 1 % model phase velocity
                  hold(hca,'on')
                  line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
                  hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
                  irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
                  hold(hca,'off')
                end
                if 1 % observed phase velocity
                  hold(hca,'on')
                  hlines = plot(hca,obs_eh_xvec,obs_velocity*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
                  irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
                  hold(hca,'off')
                end  
                colormap(hca,cn.cmap('blue_red'))
              end
              if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
                hca = h(isub); isub = isub + 1;
                units_scale_2 = 1e6;
                plot_EDI = ts_edi_flux180.data/units_scale_2;
                plot_model = ts_mod_flux_180.data/units_scale_2;
                plot(hca,x_vec,plot_model,ts_edi_flux180.time-t0,plot_EDI);                 
                hca.YLabel.String = {'flux',sprintf('(10^%g cm^{-2}s^{-1})',log10(units_scale_2))};
                legend(hca,'model','EDI','model displaced','location','eastoutside')
                text(hca,0,hca.YLim(2),'180^o','verticalalignment','top')
                hca.YLim(1) = 0;
                irf_legend(hca,{sprintf('Q_X=%.2f, Q_A=%g',Q_Flux_xcorr,Q_Flux_area)},[0.01 0.1],'color',[0 0 0])
              end
              if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0
                hca = h(isub); isub = isub + 1;
                units_scale_2 = 1e6;
                plot_EDI = ts_edi_flux0.data/units_scale_2;
                plot_model = ts_mod_flux_0.data/units_scale_2;
                plot(hca,x_vec,plot_model,ts_edi_flux0.time-t0,plot_EDI);                 
                hca.YLabel.String = {'flux',sprintf('(10^%g cm^{-2}s^{-1})',log10(units_scale_2))};
                legend(hca,'.model','EDI','model displaced','location','eastoutside')
                text(hca,0,hca.YLim(2),'0^o','verticalalignment','top')
                hca.YLim(1) = 0;
              end
              if 0 % F tot
                hca = h(isub); isub = isub + 1;
                pcolor(hca,X,V*1e-6,Ftot);
                shading(hca,'flat')
                hca.YLabel.String = 'v (km/s)';
                hcb = colorbar('peer',hca);
                hcb.YLabel.String = 'F';
                hca.CLim = hca.CLim(2)*[-1 1]; 
                colormap(hca,cn.cmap('blue_red'))
              end
              if 0 % sum FVdv, bulk velocity
                hca = h(isub); isub = isub + 1;
                plot(hca,x_vec,sumFVdv./sumFdv*1e-3);  
                hca.YLabel.String = {'v','(km/s)'};
              end          

              axes_width = 0.4;
              axes_x = 0.1;
              for ipanel = 1:npanels
                hca = h(ipanel);
                hca.Position(1) = axes_x;
                hca.Position(3) = axes_width;
                hca.Position(4) = hca.Position(4)*1.3;
                hca.Position(3);
                hca.XLim = x_vec([1 end]);
              end
              for ipanel = 1:(npanels-1)
                hca = h(ipanel);
                hca.XTickLabel = [];  
              end
              h(end).XLabel.String = sprintf('time since %s (s)',tint_phi(1).utc);
              linkaxes(h,'x')
              irf_plot_axis_align          
              h(1).YLim = [-49 49];

              if 1 % average over time, comaprison to FPI, add to right side of figure
                %h = subplot(nrows,3,[3 6 9 12]);
                h = subplot(nrows,2,[2 4 6 8]+0);
                h.Position = [0.65    0.65    0.3    0.3];
                isub = 1;
                if 1 % F, 
                  hca = h(isub); isub = isub + 1;
                  v_scale = 1e-3;
                  hlines = plot(hca,v_vec*v_scale*1e-3,mod_f_average,v_vec*v_scale*1e-3,mod_fmax_average,v_fpi*v_scale,f_fpi,'--');
                  hca.YLabel.String = {'f (s^1/m^4)'};
                  hca.XLabel.String = {'v (10^3 km/s)'};
                  hca.XLim = [-40 40];
                  fpi_utc = ts_f_fpi.time.utc;
                  str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
                    sprintf('-- fpi: %s',fpi_utc(1,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(2,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(3,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(4,12:23))};
                  %legend(hlines,str_lines)
                  irf_legend(hca,str_lines,[0.99 0.99])
                  str_info = {['T_{in}= [' sprintf('%g  ',T) '] eV'];...
                    ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
                    ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
                    sprintf('beta_{Schamel}=%g',beta);...
                    };
                  set(hca,'ColorOrder',zeros(10,3))
                  irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
                  if 1 % EDI velocities                
                    hold(hca,'on')
                    all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                               -v_edi_minus; -v_edi_plus]*[1 1];
                     if 1
                       plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
                       irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
                     end
                    hold(hca,'off')
                  end
                end
              end
              if 0 % flat skymap, to check how normalization is done
                %h = subplot(nrows,3,[3 6 9 12]);
                h = subplot(nrows,2,[2 4]+8);
                h.Position = [0.65    0.35    0.3    0.2];
                isub = 1;
                if 1 % F, 
                  hca = h(isub); isub = isub + 1;
                  [ax,hcb] = mms.plot_skymap(hca,ePDist1,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
                  %hca.YLabel.String = {'f','(s^1/m^4)'};
                  %hca.XLabel.String = {'v','(10^3 km/s)'};
                  str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
                    sprintf('-- fpi: %s',fpi_utc(1,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(2,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(3,12:23));...
                    sprintf('-- fpi: %s',fpi_utc(4,12:23))};
                  %legend(hlines,str_lines)
                  %irf_legend(hca,str_lines,[0.99 0.99])                  
                end
              end
              if 0 % flat skymap, units of flux, to check how normalization is done
                %%
                %h = subplot(nrows,3,[3 6 9 12]);
                h = subplot(nrows,2,[2 4]+14);
                h.Position = [0.65    0.05    0.3    0.2];
                isub = 1;                
                hca = h(isub); isub = isub + 1;
                if 1
                echannel = 17;
                E_fpi = mean(ePDist1.tlim(tint_phi).depend{1}(:,echannel),1);
                f_485 = squeeze(mean(ePDist1.convertto('s^3/m^6').tlim(tint_phi).data(:,echannel,:,:)));
                ang_polar = ePDist1.tlim(tint_phi).depend{3};
                ang_azim = mean(ePDist1.tlim(tint_phi).depend{2},1);
                d_azim = pi/16;
                d_pol = 1 - cosd(180/16);
                d_vel_vec = linspace(v_edi_minus,v_edi_plus,100);
                d_vel = trapz(d_vel_vec,d_vel_vec.^2); % [v^2dv] = m3/s3
                v_fpi = sqrt(2*units.e*E_fpi./units.me); % m/s
                f_vol = d_vel*d_azim*d_pol;
                flux_fpi_485 = f_485*v_fpi*f_vol; % m-2s-1, (to make into cm-2s-1, multiply with 1e-4)
                flux_scale = 1e6;
                surf(hca,ang_azim,ang_polar,flux_fpi_485'*1e-4/flux_scale); % cm-2s-1
                hcb = colorbar('peer',hca);
                hcb.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(flux_scale));
                hca.YLabel.String = 'Polar angle';
                hca.XLabel.String = 'Azimuthal angle';
                hca.Title.String = 'Where electrons are going';
                view(hca,[0 0 1])
                hca.XLim = [0 360];
                hca.YLim = [0 180];
                %[ax,hcb] = mms.plot_skymap(hca,ePDist1.convertto('s^3/m^6'),'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
                else
                  %%
                  [ax,hcb] = mms.plot_skymap(hca,ePDist1.vd3v,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
                end
                %hca.YLabel.String = {'f','(s^1/m^4)'};
                %hca.XLabel.String = {'v','(10^3 km/s)'};
                str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
                  sprintf('-- fpi: %s',fpi_utc(1,12:23));...
                  sprintf('-- fpi: %s',fpi_utc(2,12:23));...
                  sprintf('-- fpi: %s',fpi_utc(3,12:23));...
                  sprintf('-- fpi: %s',fpi_utc(4,12:23))};
              end
              if doPrint
                cn.print(sprintf('schamel_2F_vph%g_ntot%g_R%g_T1%g_T2%g_vd1%g_vd2%g_beta%g_phishift%g',vph,ntot*1e-6,R,T(1),T(2),vd(1)*1e-3,vd(2)*1e-3,beta,phi_shift))
              end
              
              if 0 % average over time, comaprison to FPI
                figure(34)
                hca = subplot(nrows,3,[3 6 9]);
                clear h;
                nrows = 1;
                ncols = 1;
                npanels = nrows*ncols;
                for ip = 1:npanels
                  h(ip) = subplot(nrows,ncols,ip);
                end
                isub = 1;
                if 1 % F, 
                  hca = h(isub); isub = isub + 1;
                  v_scale = 1e-3;
                  hlines = plot(hca,v_vec*v_scale*1e-3,mod_f_average,v_vec*v_scale*1e-3,mod_fmax_average,v_fpi*v_scale,f_fpi,'--');
                  hca.YLabel.String = {'f','(s^1/m^4)'};
                  hca.XLabel.String = {'v','(10^3 km/s)'};
                  hca.XLim = [-40 40];
                  str_lines = {'f_{mod}';'f_{mod,\phi=0}';'-- fpi';'-- fpi';'-- fpi';'-- fpi'};
                  %legend(hlines,str_lines)
                  irf_legend(hca,str_lines,[0.99 0.99])
                  str_info = {['T_{in}= [' sprintf('%g  ',T) '] eV'];...
                    ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
                    ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
                    sprintf('beta_{Schamel}=%g',beta);...
                    };
                  set(hca,'ColorOrder',zeros(10,3))
                  irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
                  if 1 % EDI velocities                
                    hold(hca,'on')
                    all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                               -v_edi_minus; -v_edi_plus]*[1 1];
                     if 1
                       plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
                       irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
                     end
                    hold(hca,'off')
                  end
                end
                if doPrint 
                cn.print(sprintf('schamel_2F_average_vph%g_ntot%g_R%g_T1%g_T2%g_vd1%g_vd2%g_beta%g',ntot*1e-6,R,T(1),T(2),vd(1)*1e-3,vd(2)*1e-3,beta))
              end
              end
            end
          end
        end
      end    
    end
  end
end

