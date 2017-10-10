function Bz = L9_time_petscheck_Bz(t)

if t>0
  Bz = (units.c/vA)*mr.L9_time_petscheck_Er(t);
else
  Bz = -(units.c/vA)*mr.L9_time_petscheck_Er(t);
end
