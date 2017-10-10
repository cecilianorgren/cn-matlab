function [mB,cmaxs,mE,mintEdt,ks,ns,Bs,B0s,cs,mva_v1,mva_v2,mva_v3,mva_l123,EB,tints,ffilts,vs] = moving_average(E,B,flh,flim,nTlh,n,v_range)
% out = tool.moving_average(irf_tlim(diE,tint),irf_tlim(diB,tint),irf_tlim(flh,tint),flim,nTlh,n,v);
% c_eval('out = tool.moving_average(irf_tlim(diE?,tint),irf_tlim(diB?,tint),irf_tlim(flh?,tint),flim,nTlh,n,v);',sc);
%
%   out = {mB,cmaxs,mE,mintEdt,ks,ns,Bs,B0s,cs,mva_v1,mva_v2,mva_v3,mva_l123,EB,tints,ffilts};
tint = [E(1,1) E(end,1)];

%f_lim = 0.8;
f_filt = [flh(:,1) flh(:,2)*flim];

if 1 % constant/average flh
    % number of wave periods to include in each moving bin, try 5?
    % just take it for an average flh.. flh = 60 gives T = nTlh*1/flh
    %nTlh = 5; % number of wave periods
    flh_av = mean(flh(:,2));
    T = nTlh/flh_av;
    %nT = diff(tint_a)/T;
    t0s = tint(1) + 0:0.5*T:tint(2); % start of each moving time interval
    nT = numel(t0s);
    if 0
        %%
        irf_plot(flh); hold on
        hhh=irf_plot([t0s' t0s'./t0s'*flh_av])
        hhh.Marker = 'x';
        hold off
    end
else % different/moving flh
    t0s(1) = tint(1);
    tt = 0;
    while 1
        tt=tt+1;
        flh0s = flh(flh(:,1)-tt0s(tt),1);
        t0s(tt) = 000;
    end    
end

angles=1:3:360;

f_filt_av = flh_av;
%B = irf_filt(B(:,1:4),f_filt_av,0);

% initialize variables
ts = nan(nT,4);
ns = nan(nT,4);
Bs = nan(nT,4);
cs = nan(nT,numel(angles));
cmaxs = nan(nT,2);
mB = nan(nT,2);
mE = nan(nT,2);
mintEdt = nan(nT,2);
B0s = nan(nT,2);
Ts = nan(nT,2);
tints = nan(nT,3);
mva_v1 = nan(nT,4);
mva_v2 = nan(nT,4);
mva_v3 = nan(nT,4);
mva_l123 = nan(nT,4);
ffilts = nan(nT,2); 
vs = nan(nT,2); % velocities

kmvas = zeros(nT,4);
mvaang = zeros(nT,2);
l12 = zeros(nT,2);
l23 = zeros(nT,2);    

for kk = 1:nT  
    flh_loc = irf_tlim(flh,t0s(kk)+0.1*[-1 1]);
    T_loc = nTlh/mean(flh_loc(:,2)); 
    tint_loc = t0s(kk) + 0.5*T_loc*[-1 1];    
    t_mean = t0s(kk);
    Ts(kk,:) = [t_mean T_loc];
    tints(kk,:) = [t_mean tint_loc];

    % take out lower hybrid frequency for filtering.
    f_filt_loc = cn_average2(irf_tlim(f_filt,tint_loc),1);
    n_loc = cn_average2(irf_tlim(n,tint_loc),1);
    locE = irf_tlim(E,tint_loc);
    locB = irf_tlim(B(:,1:4),tint_loc);
    if any([~isempty(locE) ~isempty(locB)])
        try
        [x,y,z,corr_dir,intEdt,Bz,B0,dEk,dEn,Ek,En]=irf_match_phibe_dir(locB,locE,angles,f_filt_loc);    
        mB(kk,:) = [t_mean mean(abs(Bz(:,2)))];

        cmaxs(kk,:) = [t_mean max(corr_dir(:,1))];
        i_dir = find(corr_dir(:,1)==max(corr_dir(:,1)));    
        mE(kk,:) = [t_mean mean(abs(dEk(2:end,i_dir+1)))];
        mintEdt(kk,:) = [t_mean mean(abs(intEdt(2:end,i_dir+1)))];

        ks(kk,:) = [t_mean x(i_dir,:)];
        ns(kk,:) = [t_mean y(i_dir,:)];
        Bs(kk,:) = [t_mean z(i_dir,:)];

        B0s(kk,:) = [t_mean B0];
        cs(kk,:) = corr_dir'';

        [~,l,v]=irf_minvar(irf_tlim(locE,tint_loc));
        mva_v1(kk,:) = [t_mean v(1,:)];
        mva_v2(kk,:) = [t_mean v(2,:)];
        mva_v3(kk,:) = [t_mean v(3,:)];
        mva_l123(kk,:) = [t_mean l];
        ffilts(kk,:) = [t_mean f_filt_loc];
        
        %% Velocity and density
        %n=linspace(0.1,5,100);
        v_vec=linspace(20,1000,30);
        try 
            [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n_loc,v_vec);
            i_v=find(corr_v(:)==min(corr_v(:)));
            velocity=v_vec(i_v);
        catch
            velocity = NaN;
        end        
        vs(kk,:) = [t_mean velocity];
      
        %disp([num2str(kk) '/' num2str(nT) ', Cmax = ' num2str(max(corr_dir(:,1))) ', v = ' num2str(velocity) ' km/s, ' irf_time(tint_loc(1),'epoch>utc') ' ' irf_time(tint_loc(2),'epoch>utc')])
        if max(corr_dir(:,1))>0.6%any(kk == 10:18)
            %%
            h = irf_plot(4);
            isub = 1;  
            
            hca = h(isub); isub = isub+1;
            irf_plot(hca,Bz)
            hca.YLabel.String = '\delta B_{||} [nT]';
            title(hca,['N = ' num2str(kk) ', f_{filt} = ' num2str(f_filt_loc,'%.1f') ' Hz , C = ' num2str(max(corr_dir(:,1)))])
            
            hca = h(isub); isub = isub+1;
            irf_plot(hca,intEdt(:,[1 i_dir+1]))
            hca.YLabel.String = '\int Edt [s mV/m]';            
            
            hca = h(isub); isub = isub+1;
            irf_plot(hca,[intEdt(:,1) intEdt(:,i_dir+1)/max(intEdt(:,i_dir+1))]); hold(hca,'on')
            irf_plot(hca,[Bz(:,1) Bz(:,2)/max(Bz(:,2))]); hold(hca,'off')
            hca.YLabel.String = '|\phi|/\phi_{max}';
            irf_legend(hca,{'\int E dt','B_{||}'},[0.98, 0.95]);
            irf_zoom(hca,'y',1.1*[-1 1])

                        
            hca = h(isub); isub = isub+1;
            irf_plot(hca,phi_E(:,[1 i_v+1])); hold(hca,'on')
            irf_plot(hca,phi_B); hold(hca,'off')
            hca.YLabel.String = '\phi [V]';
            irf_legend(hca,{['\phi_E (v =' num2str(velocity,'%.0f') ' km/s)'],'\phi_B'},[0.98, 0.95]);
            irf_zoom(hca,'y')
            
            pause            
        end
        catch
            disp(['some error in kk = ' num2str(kk)])
        end
    else
        disp(['kk = ' num2str(kk) ' has some empty field'])
    end
end

EB = irf_multiply(1e-3/1e-9,mE,1,mB,-1);

%scpmNe=irf_resamp(scpNe,ks(:,1));
%mv = [ks(:,1) 1e-3*B0s(:,2)*1e-9.*mB(:,2)*1e-9./(units.e*units.mu0.*scpmNe(:,2)*1e6.*mintEdt(:,2)*1e-3)];

out = {mB,cmaxs,mE,mintEdt,ks,ns,Bs,B0s,cs,mva_v1,mva_v2,mva_v3,mva_l123,EB,tints,ffilts,v};


if 0 % plot all the maximum correlation fields
        nTints = numel(tints);
    nFigures = ceil(nTints/10);
    iTint = 0;
    nPanels = 10;
    sc = 2;
    for iFigure = 1:nFigures;
        figure; 
        h = irf_plot(10);
        isub = 1;
        for iPanel = 1:nPanels
            iTint = iTint+1; if iTint>nTints; break; end
            hca = h(isub); isub = isub+1;
            %c_eval('plotE = irf_tlim(diE?,tints(iTint,2:3));',sc)
            c_eval('irf_plot(hca,irf_tlim(diE?,tints(iTint,2:3)));',sc)

        end
        pause
    end
end