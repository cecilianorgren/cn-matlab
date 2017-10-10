% Called by ExB.run and makes various plot defined by toplot.
switch toplot
    case 1 % Test plot what to expect, debug
        fig = figure(toplot);
        for k = 1:6; h(k) = subplot(3,2,k); end
        isub = 1;
        irf_colormap('poynting');

        [RS,ZS] = meshgrid(rSurf,zSurf);
        Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Etot = sqrt(Er.^2+Ez.^2);

        if 1 % Illustrate the amplitude of the radial electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m    
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
            %caxis(hca,20*[-1 1]*1e3)    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Illustrate the amplitude of the parallel electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
            %caxis(hca,20*[-1 1]*1e3)    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Illustrate the amplitude of the parallel electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Etot*1e3)    
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
            %caxis(hca,20*[-1 1]*1e3)    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Illustrate the expected ExB drift       
            ExBdrift = - Er / (B0*1e-9); % m/s 
            ExBdrift = ExBdrift*1e-3; % km/s

            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,ExBdrift)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha3 = colorbar('peer',hca); ylabel(cha3,'km/s')
            %caxis(hca,20*[-1 1]*1e3)

            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 0 % Illustrate the grid           
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,RS,ZS,RS)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha4 = colorbar('peer',hca); ylabel(cha4,'r [km]')    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 0 % Illustrate the grid           
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,RS,ZS,ZS)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha5 = colorbar('peer',hca); ylabel(cha5,'z [km]')       
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end

        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
        title(h(1),strTitle') 
    case 2 % Plot small overview
        fig = figure(toplot);
        if isInvisible; set(fig,'visible','off'); end
        for k = 1:6; h(k) = subplot(2,3,k); end
        set(fig,'position',[1 1 1680 955]); % full screen
        isub = 1;
        irf_colormap('poynting');

        [RS,ZS] = meshgrid(rSurf,zSurf);
        Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        ExBdrift = - Er / (B0*1e-9); % m/s 
        ExBdrift = ExBdrift*1e-3; % km/s

        if 1 % Illustrate the amplitude of the radial electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m  
            Elim = max(max(Er*1e3));
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
            caxis(hca,Elim*[-1 1])    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Illustrate the amplitude of the parallel electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
            %caxis(hca,20*[-1 1]*1e3)    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Illustrate the expected ExB drift           
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,ExBdrift)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)

            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
            title(hca,'ExB-drift')
        end
        if 1 % 1D starting velocities, vx vy vz
            hca = h(isub); isub = isub + 1;
            hist(hca,[vx0*1e-3,vy0*1e-3,vz0*1e-3],30);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v [10^3 km/s]'); ylabel(hca,'# of particles');
            legend(hca,'v_x','v_y','v_z')
        end  
        if 1 % Un normalized number of passes per bin
            hca = h(isub); isub = isub + 1;
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid',zGrid,surfa,numbrz')
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            cha2 = colorbar('peer',hca);
            ylabel(cha2,'Particles per bin');
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end
        if 0 % Normalized number of passes per bin
            hca = h(isub); isub = isub + 1;
            normnumbrz = repmat(rSurf,numel(zSurf),1);
            surf(hca,rGrid',zGrid,surfa,numbrz'./normnumbrz)
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            cha2 = colorbar('peer',hca);
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end
        if 1 % Plot average azimuthal velocity
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid',zGrid,surfa,vazmat'*1e-3)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'km/s');
            %vlim = max(max(abs(vazmat)));
            %caxis(hca,20*[-1 1]*1e3)
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)           
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
            title(hca,'Average azimuthal velocity')
        end
        if 0 % Starting positions
            hca = h(isub); isub = isub + 1;    
            plot3(hca,x0,y0,z0,'.'); hold(hca,'on');
            quiver3(hca,x0,y0,z0,vx0,vy0,vz0); hold(hca,'off');
            title(hca,'Starting positions')
            view(hca,[0 0 1])
            xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
            set(hca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
        end         
        if 0 % v_perp
            hca = h(isub); isub = isub + 1;
            hist(hca,sqrt(vx0.^2+vy0.^2)*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Perpendicular starting velocities')
            xlabel(hca,'v_{\perp} [10^3 km/s]'); ylabel(hca,'# of particles');
        end 
        if 0 % v_tot
            hca = h(isub); isub = isub + 1;
            hist(hca,sqrt(vx0.^2+vy0.^2+vz0.^2)*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Total starting velocities')
            xlabel(hca,'v_{tot} [10^3 km/s]'); ylabel(hca,'# of particles');
        end 
        if 0 % v_x
            hca = h(isub); isub = isub + 1;
            hist(hca,vx0*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_x [10^3 km/s]'); ylabel(hca,'# of particles');
        end    
        if 0 % v_y
            hca = h(isub); isub = isub + 1;
            hist(hca,vy0*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_y [10^3 km/s]'); ylabel(hca,'# of particles');
        end    
        if 0 % v_z
            hca = h(isub); isub = isub + 1;
            hist(hca,vz0*1e-3,100); hold(hca,'on');
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
            f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
            f = @(v,vt) (1/vt).*exp(-(v/vt).^2);
            %vt = 50e3; % km/s
            [bincounts,binpositions] = hist(hca,vz0,100); hold(hca,'on');
            binwidth = binpositions(2) - binpositions(1);
            histarea = binwidth*sum(bincounts);
            v = linspace(0,4*vtpar,100)*1e-3;
            plot(hca,v,f(v,vtpar*1e-3)*histarea*2e-5,'r'); hold(hca,'off');

            %plot(hca,v,nParticles/10*f(v,vt*1e-3));
             hold off
        end 
        if 0 % bounce statistics
            hca = h(isub); isub = isub + 1;
            maxBounce = max(nBounce);
            [nnn,~] = histc(nBounce,0:maxBounce);    
            bar(hca,0:maxBounce,nnn);
            fractionBounce = sum(nnn(2:end))/sum(nnn);
            bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];    
            title(hca,['Bouncing particles: ' bounceStr ])
            xlabel(hca,'# of bounces'); ylabel(hca,'# of particles');    
        end   
        if 0 % overshoot statistics
            hca = h(isub); isub = isub + 1;    
            maxOver = max(nOver);
            [nnn,~] = histc(nOver,0:maxOver);    
            bar(hca,0:maxOver,nnn);    
            title(hca,'Overshoot')
            xlabel(hca,'# of overshoots'); ylabel(hca,'# of particles');
        end   
        maxBounce = max(nBounce);
        [nnn,~] = histc(nBounce,0:maxBounce);    
        %bar(hca,0:maxBounce,nnn);
        fractionBounce = sum(nnn(2:end))/sum(nnn);
        bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];            
        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles) ', ' bounceStr],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
        title(h(2),strTitle')          
    case 3 % Plot xyz binning
        fig = figure(toplot);
        if isInvisible; set(fig,'visible','off'); end
        for k = 1:6
            h(k) = subplot(2,3,k);
        end
        set(fig,'position',[1 1 1400 820]); % not full screen
        isub = 1;
        irf_colormap('poynting');        

        if 0 % Illustrate the expected ExB drift           
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,ExBdrift)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1])

            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
            title(hca,'ExB-drift')
        end
        if 1 % Plot average azimuthal velocity, rz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,vazmat'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca); 
            ylabel(chca,'10^3 km/s');
            vmax = max(max(vazmat));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,'Azimuthal velocity')
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 1 % Plot average azimuthal velocity, xy-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,vazmat2xy'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xy));
            caxis(hca,vmax*[-1 1]*1e-6)            
            title(hca,['Azimuthal velocity, z = [' num2str(zGrid(1)) ' ' num2str(zGrid(end)) '] km']')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 1 % Plot average azimuthal velocity, xy-plane, center about z = 0
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,vazmat2xyCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2xyCenter)));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,['Azimuthal velocity, z = [' num2str(zGrid(zind(1))) ' ' num2str(zGrid(zind(end))) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 1 % Plot average azimuthal velocity, xz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);    
            surf(hca,xGrid,zGrid,surfa,vazmat2xz'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xz));
            caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Azimuthal velocity')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1]) 
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 1 % Plot average azimuthal velocity, xz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);    
            
            surf(hca,xGrid,zGrid,surfa,vazmat2xzCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2xzCenter)));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,['Azimuthal velocity, z = [' num2str(yGrid(yind(1))) ' ' num2str(yGrid(yind(end))) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1]) 
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 0 % Plot average azimuthal velocity, yz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,ny+1);    
            surf(hca,yGrid,zGrid,surfa,vazmat2yz'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'y')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2yz)));
            caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Azimuthal velocity')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])   
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 1 % Plot mass
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,ny+1);    
            surf(hca,yGrid,zGrid,surfa,numbxy)    
            shading(hca,'flat')
            xlabel(hca,'y')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'mass')
            %vmax = abs(max(max(vazmat2yz)));
            %caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Mass distribution')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])   
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 0 % Un normalized number of passes per bin
            hca = h(isub); isub = isub + 1;
            surf(hca,rGrid',zGrid,surfa,numbrz')
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            cha2 = colorbar('peer',hca);
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end      

        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
        title(h(2),strTitle') 
    case 4 % 
        fig = figure(toplot);
        if isInvisible; set(fig,'visible','off'); end
        for k = 1:8; h(k) = subplot(2,4,k); end
        set(fig,'position',[1 1 1680 955]); % full screen
        isub = 1;
        irf_colormap('poynting');

        [RS,ZS] = meshgrid(rSurf,zSurf);
        Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        ExBdrift = - Er / (B0*1e-9); % m/s 
        ExBdrift = ExBdrift*1e-3; % km/s

        if 0 % Illustrate the amplitude of the radial electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m  
            Elim = max(max(Er*1e3));
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
            caxis(hca,Elim*[-1 1])    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 0 % Illustrate the amplitude of the parallel electric field               
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
            %caxis(hca,20*[-1 1]*1e3)    
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 0 % Illustrate the expected ExB drift           
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,ExBdrift)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)

            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
            title(hca,'ExB-drift')
        end
        if 0 % 1D starting velocities, vx vy vz
            hca = h(isub); isub = isub + 1;
            hist(hca,[vx0*1e-3,vy0*1e-3,vz0*1e-3],30);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v [10^3 km/s]'); ylabel(hca,'# of particles');
            legend(hca,'v_x','v_y','v_z')
        end  
        if 0 % Un normalized number of passes per bin
            hca = h(isub); isub = isub + 1;
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid',zGrid,surfa,numbrz')
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            cha2 = colorbar('peer',hca);
            ylabel(cha2,'Particles per bin');
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end
        if 1 % Normalized number of passes per bin
            hca = h(isub); isub = isub + 1;
            normnumbrz = repmat(rSurf,numel(zSurf),1);
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid',zGrid,surfa,numbrz'./normnumbrz)
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            title(hca,'#passes per bin/ r')
            cha2 = colorbar('peer',hca);
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end
        if 0 % Plot average azimuthal velocity
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid',zGrid,surfa,vazmat'*1e-3)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'km/s');
            vlim = max(max(abs(ExBdrift)));
            %caxis(hca,20*[-1 1]*1e3)
            caxis(hca,vlim*[-1 1]*2)    
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
            title(hca,'Average azimuthal velocity')
        end
        if 1 % Starting positions
            hca = h(isub); isub = isub + 1;    
            plot3(hca,x0,y0,z0,'.'); hold(hca,'on');
            %pb = patch(x0,y0,z0,'b'); hold(hca,'on');
            %alpha(pb,0.1)
            quiver3(hca,x0,y0,z0,vx0,vy0,vz0); hold(hca,'off');
            title(hca,'Starting positions')
            view(hca,[0 0 1])
            xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
            set(hca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
        end   
        if 0 % Starting positions for overshoot particles
            hca = h(isub); isub = isub + 1;    
            hhh = plot3(hca,x0overshoot*1e-3,y0overshoot*1e-3,z0overshoot*1e-3,'.'); hold(hca,'on');
            %set(hhh,'FaceAlpha',0.5)
            %quiver3(hca,x0overshoot,y0overshoot,z0overshoot,vx0overshoot,vy0overshoot,vz0overshoot); hold(hca,'off');
            title(hca,'Starting positions')
            view(hca,[0 0 1])
            xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
            set(hca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
        end   
        if 0 % v_x
            hca = h(isub); isub = isub + 1;
            hist(hca,vx0*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_x [10^3 km/s]'); ylabel(hca,'# of particles');
        end    
        if 0 % v_y
            hca = h(isub); isub = isub + 1;
            hist(hca,vy0*1e-3,100);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_y [10^3 km/s]'); ylabel(hca,'# of particles');
        end    
        if 1 % v_z
            hca = h(isub); isub = isub + 1;
            hist(hca,vz0*1e-3,30); hold(hca,'on');
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting velocities')
            xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
            f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
            f = @(v,vt) (1/vt/sqrt(pi)).*exp(-((v+veh*1e-3)/vt).^2);
            %vt = 50e3; % km/s
            [bincounts,binpositions] = hist(hca,vz0,30); hold(hca,'on');
            binwidth = binpositions(2) - binpositions(1);
            histarea = binwidth*sum(bincounts);
            v = linspace(-2,2,100)*vtpar*1e-3-veh*1e-3;
            plot(hca,v,f(v,vtpar*1e-3)*histarea*1e-3,'r'); hold(hca,'off');

            %plot(hca,v,nParticles/10*f(v,vt*1e-3));
             hold off
        end
        if 1 % v_perp
            hca = h(isub); isub = isub + 1;
            hist(hca,sqrt(vx0.^2+vy0.^2)*1e-3,30);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Perpendicular starting velocities')
            xlabel(hca,'v_{\perp} [10^3 km/s]'); ylabel(hca,'# of particles');
        end 
        if 1 % v_tot
            hca = h(isub); isub = isub + 1;
            hist(hca,sqrt(vx0.^2+vy0.^2+vz0.^2)*1e-3,30);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Total starting velocities')
            xlabel(hca,'v_{tot} [10^3 km/s]'); ylabel(hca,'# of particles');
        end 
        if 1 % bounce statistics
            hca = h(isub); isub = isub + 1;
            maxBounce = max(nBounce);
            [nnn,~] = histc(nBounce,0:maxBounce);    
            bar(hca,0:maxBounce,nnn);
            fractionBounce = sum(nnn(2:end))/sum(nnn);
            bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];    
            title(hca,['Bouncing particles: ' bounceStr ])
            xlabel(hca,'# of bounces'); ylabel(hca,'# of particles');    
        end   
        if 0 % overshoot statistics
            hca = h(isub); isub = isub + 1;    
            maxOver = max(nOver);
            [nnn,~] = histc(nOver,0:maxOver);  
            nOvershootingParticles = sum(nnn(2:end));
            nTotOvershoots = sum(nnn.*tocolumn(0:(numel(nnn)-1)));
            bar(hca,0:maxOver,nnn);    
            titleStr = {['Overshoot, sum(all but 0 overshoot) = ' num2str(nOvershootingParticles)],...
                         ['Total number of overshoots = ' num2str(nTotOvershoots)]};
            title(hca,titleStr)
            xlabel(hca,'# of overshoots'); ylabel(hca,'# of particles');
        end  
         if 1 % v_z
            hca = h(isub); isub = isub + 1;
            hist(hca,vz0flow*1e-6,30); hold(hca,'on');
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Starting overshoot velocities')
            xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
            %f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
            f = @(v,vt) (abs(v)/vt/sqrt(pi)).*exp(-((v+veh*1e-3)/vt).^2);
            %vt = 50e3; % km/s
            [bincounts,binpositions] = hist(hca,vz0flow,30); hold(hca,'on');
            binwidth = binpositions(2) - binpositions(1);
            histarea = binwidth*sum(bincounts);
            v = linspace(-3,3,100)*vtpar*1e-3-veh*1e-3;
            % something wrong with normalization for this 'fit'
            %plot(hca,v,f(v,vtpar*1e-3)*histarea*1e-8,'r'); hold(hca,'off');
            
            %plot(hca,v,nParticles/10*f(v,vt*1e-3));
             hold off
        end
        if 1 % v_perp
            hca = h(isub); isub = isub + 1;
            hist(hca,sqrt(vx0(nParticlesBox+1:end).^2+vy0(nParticlesBox+1:end).^2)*1e-6,30);
            set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
            title(hca,'Perpendicular overshoot velocities')
            xlabel(hca,'v_{\perp} [10^3 km/s]'); ylabel(hca,'# of particles');
        end 
        maxBounce = max(nBounce);
        [nnn,~] = histc(nBounce,0:maxBounce);    
        %bar(hca,0:maxBounce,nnn);
        fractionBounce = sum(nnn(2:end))/sum(nnn);
        bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];            
        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles) ', ' bounceStr],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
        title(h(2),strTitle')  
    case 5 % Plot current
        fig = figure(toplot);
        if isInvisible; set(fig,'visible','off'); end
        for k = 1:6
            h(k) = subplot(2,3,k);
        end
        set(fig,'position',[1 1 1400 820]); % not full screen
        isub = 1;
        irf_colormap('poynting');        

        if 1 % Sum(m*v_az), rz
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid',zGrid,surfa,sumMVrz'*1e-3)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'km/s');
            vlim = max(max(abs(sumMVrz)));
            %caxis(hca,20*[-1 1]*1e3)
            caxis(hca,vlim*[-1 1]*1e-3)    
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
            title(hca,'Sum(m*v_az)')
        end        
        if 1 % Sum(m), rz
            hca = h(isub); isub = isub + 1;
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid',zGrid,surfa,numbrz')
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')    
            cha2 = colorbar('peer',hca);
            ylabel(cha2,'mass');
            title(hca,'sum(m)')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
        end
        if 1 % Sum(m*v_az)/sum(m)
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid',zGrid,surfa,vazmat'*1e-3)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'km/s');
            vlim = max(max(abs(vazmat*1e-3)));
            %caxis(hca,20*[-1 1]*1e3)
            caxis(hca,vlim*[-1 1])    
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
            title(hca,'sum(m*v_ax)/sum(m)')
        end  
        if 1 % Plot mass
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,ny+1);    
            surf(hca,yGrid,zGrid,surfa,numbyz')    
            shading(hca,'flat')
            xlabel(hca,'y')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'mass')
            %vmax = abs(max(max(vazmat2yz)));
            %caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Mass distribution')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])   
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end 
        if 1 % Plot average azimuthal velocity, xy-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,numbxyCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            %vmax = max(max(vazmat2xy));
            %caxis(hca,vmax*[-1 1]*1e-6)            
            title(hca,['Mass, z = [' num2str(zGrid(1)) ' ' num2str(zGrid(end)) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 1 % Mass distribution
            hca = h(isub); isub = isub + 1;
            hist(hca,m,50)                                    
            title(hca,'Mass distribution')
            xlabel(hca,'mass'); ylabel(hca,'# of particles');    
        end   
        if 0 % Plot average azimuthal velocity, xy-plane, center about z = 0
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,vazmat2xyCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xyCenter));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,['Azimuthal velocity, z = [' num2str(zGrid(zind(1))) ' ' num2str(zGrid(zind(end))) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 0 % Plot average azimuthal velocity, xz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);    
            surf(hca,xGrid,zGrid,surfa,vazmat2xz'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xz));
            caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Azimuthal velocity')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1]) 
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
   

        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
       % title(h(2),strTitle')     
    case 6 % Plot xyz binning
        [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
        Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        Er = sqrt(Ex.^2+Ey.^2);
        ExB_y = - Ex/(B0*1e-9); % m/s 
        ExB_x = + Ey/(B0*1e-9); % m/s 
        ExB_y = ExB_y*1e-3; % km/s
        ExB_x = ExB_x*1e-3; % km/s
        
        [RS,ZRS] = meshgrid(rSurf,zSurf);
        Ezr = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZRS/lz).^2)*1e-3; % V
        Erz = ZRS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZRS/lz).^2)*1e-3; % V        
        ExBdrift = - Ezr / (B0*1e-9); % m/s 
        ExBdrift = ExBdrift*1e-3; % km/s
        
        fig = figure(toplot);
        if isInvisible; set(fig,'visible','off'); end
        for k = 1:10
            h(k) = subplot(2,5,k);
        end
        set(fig,'position',[1 1 1400 820]); % not full screen
        isub = 1;
        irf_colormap('poynting');        

        if 1 % Illustrate the expected az-ExB drift, rz-plane          
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,ExBdrift)
            shading(hca,'flat')
            xlabel(hca,'r'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
            title(hca,'ExB-drift')
        end
        if 1 % Illustrate the expected y-ExB drift, xz-plane 
            meanind = fix(size(ExB_y,2)/2);
            m_ExB_y = squeeze(mean(ExB_y(meanind+(-1:1),:,:),1));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_ExB_y')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'ExB_y-drift')
        end
        if 1 % Illustrate the expected x-ExB drift, xz-plane 
            meanind = fix(size(ExB_x,2)/2);
            m_ExB_x = squeeze(mean(ExB_x(meanind+(-1:1),:,:),1));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_ExB_x')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'ExB_x-drift')
        end  
        if 1 % Illustrate the expected x-ExB drift, xy-plane 
            meanind = fix(size(ExB_x,2)/2);
            m_ExB_x = squeeze(mean(ExB_y(:,:,meanind+(-1:1)),3));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);
            surf(hca,xGrid,yGrid,surfa,m_ExB_x')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'y')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'ExB_x-drift')
        end          
        if 0 % Illustrate the expected x-ExB drift, xz-plane 
            meanind = fix(size(ExB_x,2)/2);
            m_ExB_x = squeeze(mean(ExB_x(:,meanind+(-1:1),:),1));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_ExB_x')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'km/s')
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'ExB_x-drift')
        end       
        if 0 % Plot the radial electric field in xz-plane
            Err = sqrt(Ex.^2+Ey.^2);
            meanind = fix(size(Err,2)/2);
            m_Err = squeeze(mean(Err(:,meanind+(-1:1),:),2));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_Err')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'V')
            %vlim = max(max(abs(ExBdrift)));
            %caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'E_r')
        end                
        if 0 % Plot the radial electric field in rz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);
            surf(hca,rGrid,zGrid,surfa,Ezr)
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'V')
            vlim = max(max(abs(ExBdrift)));
            %caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[0 1],'ylim',[-zlim zlim]) 
            title(hca,'E_r')
        end 
        if 0 % Plot the x electric field in xz-plane            
            meanind = fix(size(Ex,2)/2);
            m_Ex = squeeze(mean(Ex(meanind+(-1:1),:,:),1));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_Ex')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'V')
            %vlim = max(max(abs(ExBdrift)));
            %caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'E_x')
        end         
        if 0 % Plot the y electric field in xz-plane            
            meanind = fix(size(Ey,2)/2);
            m_Ey = squeeze(mean(Ey(meanind+(-1:1),:,:),1));
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);
            surf(hca,xGrid,zGrid,surfa,m_Ey')
            shading(hca,'flat')
            xlabel(hca,'x'); ylabel(hca,'z')
            cha = colorbar('peer',hca); ylabel(cha,'V')
            %vlim = max(max(abs(ExBdrift)));
            %caxis(hca,vlim*[-1 1]*2)
            axis(hca,'equal')
            grid(hca,'off'); view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
            title(hca,'E_y')
        end            
        if 1 % Plot average azimuthal velocity, rz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nr+1);    
            surf(hca,rGrid,zGrid,surfa,vazmat'*1e-3)    
            shading(hca,'flat')
            xlabel(hca,'r')
            ylabel(hca,'z')
            chca = colorbar('peer',hca); 
            ylabel(chca,'km/s');
            vlim = max(max(abs(ExBdrift)));
            caxis(hca,vlim*[-1 1]*2);
            title(hca,'Azimuthal velocity')
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
        end
        if 0 % Plot average azimuthal velocity, xy-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,vazmat2xy'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xy));
            caxis(hca,vmax*[-1 1]*1e-6)            
            title(hca,['Azimuthal velocity, z = [' num2str(zGrid(1)) ' ' num2str(zGrid(end)) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 0 % Plot total electric field in arrows
            [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
            hca = h(isub); isub = isub + 1;        
            tpx = 1:10:numel(xSurf);
            tpy = 1:10:numel(ySurf);
            tpz = 1:10:numel(zSurf);
            quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                   Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz))
            xlabel(hca,'x')
            ylabel(hca,'y')
            zlabel(hca,'z')
            title(hca,'E field')
            axis(hca,'equal')
            grid(hca,'off')            
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
        end
        if 0 % Plot x-electric field in arrows
            hca = h(isub); isub = isub + 1;        
            tpx = 1:10:numel(xSurf);
            tpy = 1:10:numel(ySurf);
            tpz = 1:10:numel(zSurf);
            quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                   Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz)*0,Ez(tpx,tpy,tpz)*0)
            xlabel(hca,'x')
            ylabel(hca,'y')
            zlabel(hca,'z')
            title(hca,'E_x field')
            axis(hca,'equal')
            grid(hca,'off')            
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
        end        
        if 0 % Plot average azimuthal velocity, xy-plane, center about z = 0
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(ny+1,nx+1);    
            surf(hca,xGrid,yGrid,surfa,vazmat2xyCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'y')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2xyCenter)));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,['Azimuthal velocity, z = [' num2str(zGrid(zind(1))) ' ' num2str(zGrid(zind(end))) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])  
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
        end
        if 0 % Plot average azimuthal velocity, xz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);    
            surf(hca,xGrid,zGrid,surfa,vazmat2xz'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = max(max(vazmat2xz));
            caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Azimuthal velocity')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1]) 
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 0 % Plot average azimuthal velocity, xz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,nx+1);    
            surf(hca,xGrid,zGrid,surfa,vazmat2xzCenter'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'x')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2xzCenter)));
            if ~isnan(vmax); caxis(hca,vmax*[-1 1]*1e-6); end
            title(hca,['Azimuthal velocity, z = [' num2str(yGrid(yind(1))) ' ' num2str(yGrid(yind(end))) '] km'])
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1]) 
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end
        if 0 % Plot average azimuthal velocity, yz-plane
            hca = h(isub); isub = isub + 1;        
            surfa = zeros(nz+1,ny+1);    
            surf(hca,yGrid,zGrid,surfa,vazmat2yz'*1e-6)    
            shading(hca,'flat')
            xlabel(hca,'y')
            ylabel(hca,'z')
            chca = colorbar('peer',hca);
            ylabel(chca,'10^3 km/s')
            vmax = abs(max(max(vazmat2yz)));
            caxis(hca,vmax*[-1 1]*1e-6)
            title(hca,'Azimuthal velocity')
            %caxis(hca,20*[-1 1]*1e3)
            %caxis(hca,vlim*[-1 1])   
            axis(hca,'equal')
            grid(hca,'off')
            view(hca,[0 0 1])
            set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
        end    

        strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
                    ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
                    };
        %title(h(2),strTitle')    
end