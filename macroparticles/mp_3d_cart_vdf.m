function out = mp_3d_cart_vdf(MP,vx_edges,vy_edges,vz_edges)

% residue/root-mean-scquare difference/whatever quality
dvx = vx_edges(2)-vx_edges(1);
dvy = vy_edges(2)-vy_edges(1);
dvz = vz_edges(2)-vz_edges(1);

[count edges mid loc] = histcn([MP.vx, MP.vy, MP.vz],vx_edges,vy_edges,vz_edges,'AccumData',MP.dn);
dvx_cms = dvx*1e5; % cm/s
dvy_cms = dvy*1e5;
dvz_cms = dvz*1e5;

dn_km3 = 1e5^3*count; % cm^-3 -> km^-3
fcart = dn_km3/(dvx*dvy*dvz);
Fobs.f = fcart;
Fobs.edges = edges;
Fobs.mid = mid;
Fobs.dv = dvx*dvy*dvz;
out = Fobs;
%out.f_units = 'typically s^3/cm^6';
out.f_units = 's^3/km^6';