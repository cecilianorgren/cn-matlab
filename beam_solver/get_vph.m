[k,w] = ginput(1);
vph = w/k*normal_frequency*normal_length;
disp(['vph = ' num2str(vph/1000) ' km/s'])