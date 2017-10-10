units=irf_units;

dedw = @(me,mp) 2*( (me/mp/16)^(1/3) ...
      - (mp/me)^(1/2)/(1-(me/mp)^(2/3)*16^(-1/3)) );
dedw(units.me,units.mp)

wk = @(me,mp) (mp/me/16)^(1/3);

wk(units.me,units.mp)

dedw(units.me,units.mp)*wk(units.me,units.mp)

Wp = @(me,mp) 16; 
We = @(me,mp) (1-(me/16/mp)^(1/3))^(-3);

Wp(units.me,units.mp)
We(units.me,units.mp)