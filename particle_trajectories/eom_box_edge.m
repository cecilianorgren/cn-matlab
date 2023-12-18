function [value, isterminal, direction] = eom_box_edge(t, x_vect,comp,boxedge)
  % integration is terminated when value changes sign
  % for this setup, value is initially negative
  value      = (x_vect(comp)-boxedge(1))*(x_vect(comp)-boxedge(2));        
  isterminal = 1;   % Stop the integration
  direction  = 0;
end