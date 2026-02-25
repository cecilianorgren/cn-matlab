for k = 1:5
    h(k) = subplot(1,5,k);
end
isub = 1;

if 1
    hca = h(isub); isub = isub + 1;
    plot3(hca,x*1e-3,y*1e-3,z*1e-3); 
    %plot3(hca,x(binz(iz))*1e-3,y(binz(iz))*1e-3,z(binz(iz))*1e-3,'g*'); hold(hca,'off');
    title(hca,'Electron trajectory')
    xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
end
if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,z*1e-3,x*1e-3,z*1e-3,y*1e-3); 
    %plot3(hca,x(binz(iz))*1e-3,y(binz(iz))*1e-3,z(binz(iz))*1e-3,'g*'); hold(hca,'off');
    title(hca,'Electron trajectory')
    xlabel(hca,'z [km]'); ylabel(hca,'x, y [km]'); zlabel(hca,'z [km]')
end
if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,z*1e-3,x*1e-3,z*1e-3,y*1e-3); 
    %plot3(hca,x(binz(iz))*1e-3,y(binz(iz))*1e-3,z(binz(iz))*1e-3,'g*'); hold(hca,'off');
    title(hca,'Electron trajectory')
    xlabel(hca,'z [km]'); ylabel(hca,'x, y [km]'); zlabel(hca,'z [km]')
end
if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,z*1e-3,th,'.')    
    title(hca,'Angular position')
    xlabel(hca,'z [km]'); ylabel(hca,'v_{\phi} [km/s]');
end
if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,z*1e-3,vaztmp*1e-6)    
    title(hca,'Azimuthal velocity')
    xlabel(hca,'z [km]'); ylabel(hca,'v_{\phi} [10^3 km/s]');
end

