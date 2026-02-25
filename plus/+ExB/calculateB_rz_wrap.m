%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^2)
modelstorun = [1 2];
toLoad = {'/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/20140324T224902_2007-08-31_1726898.mat',...
    %'/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/20140326T033835_2007-08-31_5182101.mat',...
          '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/20140324T192651_Tao2011_.mat'};
density = [0.07 0.02];
ncalc = numel(modelstorun);
vecSimBx = cell(1,ncalc);
vecSimBy = cell(1,ncalc);
vecSimBz = cell(1,ncalc);
vecSimBr = cell(1,ncalc);
vecSimfX = cell(1,ncalc);
vecSimfY = cell(1,ncalc);
vecSimfZ = cell(1,ncalc);
vecSimfR = cell(1,ncalc);
vecModBx = cell(1,ncalc);
vecModBy = cell(1,ncalc);
vecModBz = cell(1,ncalc);
vecModBr = cell(1,ncalc);
vecModfX = cell(1,ncalc);
vecModfY = cell(1,ncalc);
vecModfZ = cell(1,ncalc);
vecModfR = cell(1,ncalc);
for mm = 1:2
    modelnr = modelstorun(mm); ExB.model;
    load(toLoad{mm})
    n = density(mm);
    rbin = 9; thbin = 7; zbin = 7;
    ExB.calculateB_rz;

    vecSimBx{mm} = simBx;
    vecSimBy{mm} = simBy;
    vecSimBz{mm} = simBz;
    vecSimBr{mm} = simBr;
    vecSimfX{mm} = simRx;
    vecSimfY{mm} = simRy;
    vecSimfZ{mm} = simRz;
    vecSimfR{mm} = simRr;
    vecModBx{mm} = modBx;
    vecModBy{mm} = modBy;
    vecModBz{mm} = modBz;
    vecModBr{mm} = modBr;
    vecModfX{mm} = modRx;
    vecModfY{mm} = modRy;
    vecModfZ{mm} = modRz;
    vecModfR{mm} = modRr;
    vecRG{mm} = rg;
    vecTHG{mm} = thg;
    vecZG{mm} = zg;

end
ExB.artfig2;
if 0 % save
    
    
end