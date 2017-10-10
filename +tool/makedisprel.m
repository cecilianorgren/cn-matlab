em_tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 14.0])';
em_tint = toepoch([2008 04 22 18 06 05;2008 04 22 18 06 07])';


f_center = [0.01:0.1:1.5]';
fint = 0.1*[-1 1];
fn = numel(f_center);
csys = 'gsm';

%tint = em_tint;
sc = 3; 

save_c = zeros(fn,1);
save_v = zeros(fn,3);
save_cv = zeros(fn,1);
save_ExB = zeros(fn,3);
save_EoverB = zeros(fn,1);

for ii = 1:fn
    flim = f_center(ii)+fint;    
    if flim(1)<0.01; flim(1) = 0.02; end    
    tool.single_event
    save_phiE{ii} = phi_E(:,[1 1+i_v]);
    save_phiB{ii} = phi_B;
    save_phiEB{ii} = [phi_E(:,[1 1+i_v]) phi_B(:,2)];
    save_dEk{ii} = dEk(:,[1 1+i_dir]);
    save_dEn{ii} = dEn(:,[1 1+i_dir]);
    save_dBz{ii} = Bz;
    save_c(ii,1) = corr_dir(i_dir);
    save_v(ii,:) = velocity*x(i_dir,:);
    save_cv(ii,:) = corr_v(i_v);
    save_ExB(ii,:) = vExB;    
    save_Ekmax(ii,1) = max(save_dEk{ii}(:,2));
    save_Enmax(ii,1) = max(save_dEn{ii}(:,2));
    save_Bmax(ii,1) = max(save_dBz{ii}(:,2));    
end
save_EkoverB = save_Ekmax*1e-3./save_Bmax/1e-9;
save_EnoverB = save_Enmax*1e-3./save_Bmax/1e-9;
%
save_vabs =irf_abs(save_v,1);
save_vabs_ionframe =irf_abs(save_v-save_ExB,1);
save_ExBparv = irf_norm(save_v)*save_ExB(1,:)';
save_lam = save_vabs./(f_center*flh_loc); % km
save_lam_ionframe = save_lam; % km
save_k = 2*pi./save_lam; % 1/km
save_k_ionframe = save_k; % 1/km
save_f_ionframe = f_center - save_ExBparv./save_lam/flh_loc;
save_f_ionframe = save_vabs_ionframe.*save_k/2/pi;
c_eval('re=irf_resamp(re3,(tint),''mean'');',sc); re=re(2);
%f_ionframe = f_center - irf_abs(save_ExB,1).*save_k/flh_loc/2/pi;
if 0
%%

tint_str= [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
ax=subplot(3,1,1);
plot(f_center,save_vabs,f_center,irf_abs(save_ExB,1),f_center,save_vabs_ionframe,f_center,save_ExBparv);
ax(1).XLabel.String = '\omega/\omega_{LH}';
ax(1).YLabel.String = 'v [km/s]';
title(ax(1),tint_str)
irf_legend(ax,{'v_{ph}','v_{ExB}','v_{ph}-v_{ExB}','k\cdot v_{ExB}'},[0.02 0.95])


ax= subplot(3,1,3);
clim = 0.7;
plot(ax,save_k(save_c>clim)*re,f_center(save_c>clim),'*',save_k*re,f_center)%,...
       %save_k(save_c>clim)*re,save_f_ionframe(save_c>clim),'*',save_k*re,save_f_ionframe);
     
legend(ax(1),['C_k>' num2str(clim)])
ax(1).YLabel.String = '\omega/\omega_{LH}';
ax(1).XLabel.String = 'k\rho_e';
if ax(1).XLim(2)<2; ax(1).XLim(2)=2; end

hold(ax,'off');

subplot(3,1,2);
[ax,p1,p2] = plotyy(f_center,save_cv,f_center,save_c,'semilogy','plot');
ax(1).XLabel.String = '\omega/\omega_{LH}';
ax(1).YLabel.String = 'C_v';
ax(2).YLabel.String = 'C_k';
ax(2).YLim = [0 1];

%% Make spectrum of given time interval
specE = irf_powerfft(phi_E(:,[1 i_dir+1]),512,450,0);
specB = irf_powerfft(phi_B,512,450,0);
plot(1:900,specE,1:900,specB)
%%
specE = fft(phi_E(:,[1 i_dir+1]),1024);
plot(specE)
%%
subplot(1,2,1)
plot(f_points,irf_abs(save_v,1))
subplot(1,2,2)
plot(f_points*flh_loc,save_c)
end
    