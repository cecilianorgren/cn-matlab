for k = 1:6
    h(k) = subplot(2,3,k);
end
isub = 1;

hca = h(isub); isub = isub + 1;
plot3(hca,x*1e-3,y*1e-3,z*1e-3); hold(hca,'on');
%plot3(hca,x(binz(iz))*1e-3,y(binz(iz))*1e-3,z(binz(iz))*1e-3,'g*'); hold(hca,'off');
title(hca,'Electron trajectory')
xlabel(hca,'x'); ylabel(hca,'y'); zlabel(hca,'z')

hca = h(isub); isub = isub + 1;
plot(hca,z*1e-3,binz,z*1e-3,binr)
%title(hca,'Bin number')
xlabel(hca,'z [km]'); ylabel(hca,'Bin number');
legend(hca,'z','r')

hca = h(isub); isub = isub + 1;
plot(hca,1:numel(nrr),nrr,1:numel(nzz),nzz)
title(hca,'Passes per bin')
xlabel(hca,'Bin number (index)');
legend(hca,'r','z')

if 0
    hca = h(isub); isub = isub + 1;
    plot(hca,1:size(arz,1),arz(:,1),1:size(arz,1),arz(:,2))
    title(hca,'arz')
    %xlabel('x'); xlabel('y'); xlabel('z')
end

if 1
    hca = h(isub); isub = isub + 1;
    for k = 1:nParticles
        plot3(hca,x0(k),y0(k),z0(k),'*'); hold(hca,'on')
    end
    title(hca,'Starting positions')
    xlabel('x'); ylabel('y'); zlabel('z')
end
    

hca = h(isub); isub = isub + 1;
plot(hca,z*1e-3,vaztmp*1e-3)
title(hca,'Azimuthal velocity')
xlabel(hca,'z [km]'); ylabel(hca,'v_{\phi} [km/s]');

hca = h(isub); isub = isub + 1;
plot(hca,z*1e-3,crz)
title(hca,'crz')
xlabel('z [km]'); xlabel('y'); xlabel('z')

