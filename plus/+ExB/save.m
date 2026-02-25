% save run
% get computer id
[~, hostname] = system('hostname');
if strcmp(hostname(1:4),'spis')
    saveIn = '/home/cecilia/Research/ElectronHoleSimulation/SavedData/'; % on spis
    %disp('spis')
else    
    saveIn = '/Users/Cecilia/Research/EH/TestParticleSimulation/'; % on local
    %disp('local')
end


try
saveFilepath = [saveIn datestr(now,'yyyymmddTHHMMSS') '_' name '_' num2str(nParticles)];
save(saveFilepath,...
     'name','nBounce','x0','y0','z0','vx0','vy0','vz0',...
     'vazmat','numbrz','vazmat2yz','vazmat2xz','vazmat2xy',...
     'numbxy','numbxz','numbyz','tTot','xGrid','yGrid','zGrid','rGrid',...
     'xSurf','ySurf','zSurf','rSurf','zlim','rlim','L','lr','lz','phi0',...
     'B0','Tpar','Tper','veh','isInvisible','nx','ny','nz','nr','nParticles',...
     'vtpar','vtper','vazmat2xyCenter','vazmat2xzCenter','vazmat2yzCenter',...
     'zind','yind','xind','sumMVrz','sumMrz','sumMVxyz','sumMxyz','tend',...
     'nParticlesBox','nParticlesFlow','nParticles','z0flow','z0box',...
     'vz0box','vz0flow','m','turns','xc','yc',...
     'sumMVXxyz','sumMVYxyz','sumMVZxyz','tint','nx','ny','nz')
disp(['Saved ' saveFilepath])
catch me    
end

if 0
try
    save_run.name = name;
    save_run.nOver = nOver;
    save_run.nBounce = nBounce;
    save_run.x0 = x0;
    save_run.y0 = y0;
    save_run.z0 = z0;
    save_run.vx0 = vx0;
    save_run.vy0 = vy0;
    save_run.vz0 = vz0;
    save_run.vaz1 = vaz1;
    save_run.vazmat = vazmat;
    save_run.numbrz = numbrz;
    save('run_2014-03-01_b','save_run')
catch me
    continue
end
end
