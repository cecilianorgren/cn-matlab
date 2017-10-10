Oxygen = struct('m',16,'n',1,'t',10,'a',5,'vd',20,'d',1,'b',0);
Electrons = struct('m',0,'n',20*1e6,'t',5*1e-3,'a',1,'vd',0,'d',1,'b',0);
  [h,f,E] = whamp.plot_f(gca,Electrons,'km/s','PSDvsE','pitchangles',[0]);
   