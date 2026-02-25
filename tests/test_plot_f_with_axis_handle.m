ax=axes;
w_n=  [0.055]*1e6;       % density m^3
w_t=  [1.6];            % temperature [keV]
w_vd= [0];         % V_drift/V_term
w_a1= [1];           % T_perp/T_par
w_m=  [0];            % particle mass, 0-electrons, 1-protons, 16-oxygen
w_d=  [1];            % loss cone parameter, 1.0 = no loss cone)
w_a2= [0];           % ??? (use 0)
w_pa= [0];    % pitch angeles
plotoption=[1];
titleoption=[0];
cmp=[1];

art3.plot_f(ax,w_n(cmp),w_m(cmp),w_t(cmp),w_vd(cmp),...
                            w_d(cmp),w_a1(cmp),w_a2(cmp),w_pa,...
                            plotoption,titleoption)