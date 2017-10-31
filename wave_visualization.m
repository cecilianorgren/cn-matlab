% 'E' - wave polarization  [Ex Ey Ez] or [Ex Ez] or [Ex]
% 'kx' - X component of wave vector (default 2*pi, 1 wave period)
% 'kz' - Z component of wave vector (default 0)
% 'qm' - vector of charge/mass for each species (default [1 -1836])
% 'N' - number of particles in each species
% 'f' - wave frequency in Hz (default 1, 20 frames per s)
% 'T' - length of simulation
% 'X' - X length (default is 1)

[A,im,map]=irf_plasma_wave_visualization('E',[1 0 0.1],'k',[1 0 0.1],'f',5);