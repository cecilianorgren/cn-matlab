%% Implied electron temperature based on measured electron hole velocity and ion drift speed.
units = irf_units;
vd = @(veh,vi) (veh-vi)*(16*units.mp/units.me)^(1/3);
vte = @(veh,vi) vd(veh,vi)/0.9;
Te = @(veh,vi) cn_eV2v(vte(veh,vi),'v') ;

vi=0:440;
veh=440;
plot(vi,Te(veh,vi),200,Te(veh,200),'r*',...
     cn_eV2v(7000/0.9,'v'),130,'b*',cn_eV2v(7000,'v')/0.9,130,'g*')
xlabel('v_i')
ylabel('T_e')