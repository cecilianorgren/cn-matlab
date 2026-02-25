% save particle position and velocity to textfile
dirpath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Run1/';
if ~exist('dirpath','dir'); mkdir(dirpath); end
filenum = 1;
filepath = [dirpath 'file' num2str(filenum)]
fid = fopen(filepath,'w');
count = fprintf(fid,'%f %f %f %f %f %f/n',[x y z vx vy vz]);

% fopen(filepath) % tog j??tel?ng tid