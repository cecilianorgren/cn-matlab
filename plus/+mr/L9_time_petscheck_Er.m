function Er = L9_time_petscheck_Er(t)

Er = sin(pi*t);
Er(t<=0) = 0;
Er(t>=1) = 0;

