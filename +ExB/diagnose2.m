for k = 1:6
    h(k) = subplot(2,3,k);
end
isub = 1;

if 0
    hca = h(isub); isub = isub + 1;
    plot3(hca,x*1e-3,y*1e-3,z*1e-3); 
    %plot3(hca,x(binz(iz))*1e-3,y(binz(iz))*1e-3,z(binz(iz))*1e-3,'g*'); hold(hca,'off');
    title(hca,'Electron trajectory')
    xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
end
if 1
    hca = h(isub); isub = isub + 1;    
    plot3(hca,x0,y0,z0,'*');    
    title(hca,'Starting positions')
    view(hca,[0 0 1])
    xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
end  
if 1
    hca = h(isub); isub = isub + 1;    
    plot(hca,sqrt(x0.^2+y0.^2),z0,'*');    
    title(hca,'Starting positions')
    xlabel(hca,'r [km]'); ylabel(hca,'z [km]'); 
end  
if 1
    hca = h(isub); isub = isub + 1;
    hist(hca,vx0*1e-3,50);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_x [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 1
    hca = h(isub); isub = isub + 1;
    hist(hca,vy0*1e-3,50);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_y [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 1
    hca = h(isub); isub = isub + 1;
    hist(hca,vz0*1e-3,50); hold(hca,'on');
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
    f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
    %vt = 50e3; % km/s
    v = linspace(0,4*vtpar,100)*1e-3;
    %plot(hca,v,nParticles/10*f(v,vt*1e-3));
     hold off
end  
if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,z*1e-3,vaztmp*1e-3)    
    title(hca,'Azimuthal velocity')
    xlabel(hca,'z [km]'); ylabel(hca,'v_{\phi} [km/s]');
end

