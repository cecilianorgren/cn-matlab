function [f0,params] = get_f0(inp)

switch inp 
  case 1
    f0 = 2;
    switch f0
      case 1
        ntot = 0.04*1e6;
        R = 0.60; 
        n1 = ntot*R; 
        n2 = ntot*(1-R);
        T1 = 150;
        T2 = 2000;
        vd1 = -9000*1e3;
        vd2 = 4000*1e3; 
        
        n = [n1 n2];
        vd = [vd1 vd2];
        T = [T1 T2];
      case 2
        ntot = 0.04*1e6;
        R = 0.65; 
        n1 = ntot*R; 
        n2 = ntot*(1-R);
        T1 = 150;
        T2 = 1500;
        vd1 = -8000*1e3;
        vd2 = 6000*1e3; 
        
        n = [n1 n2];
        vd = [vd1 vd2];
        T = [T1 T2];
    end
  case 2
    ntot = 0.035*1e6;
    R = [0.4 0.25 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 3
    ntot = 0.04*1e6;
    R = 0.65; 
    n1 = ntot*R; 
    n2 = ntot*(1-R);
    T1 = 100;
    T2 = 1500;
    vd1 = -11000*1e3;
    vd2 = 6000*1e3;  
    
    n = [n1 n2];
    vd = [vd1 vd2];
    T = [T1 T2];
  case 4
    ntot = 0.04*1e6;
    R = [0.65/2 0.65/2 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 50;
    T3 = 1500;
    vd1 = -14000*1e3;
    vd2 = -6000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 5
    ntot = 0.04*1e6;
    R = [0.65/2 0.65/2 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 50;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 6
    ntot = 0.04*1e6;
    R = [0.4 0.25 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 7
    ntot = 0.04*1e6;
    R = [0.25 0.4 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = 0000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 8
    ntot = 0.035*1e6;
    R = [0.4 0.25 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 9
    ntot = 0.035*1e6;
    R = [0.3 0.35 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 10
    ntot = 0.04*1e6;
    R = [0.3 0.35 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -8000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 11
    ntot = 0.04*1e6;
    R = [0.3 0.35 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 200;
    T3 = 1500;
    vd1 = -12000*1e3;
    vd2 = -6000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 12 
    ntot = 0.04*1e6;
    R = [0.25 0.4 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 200;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -6000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 13
    ntot = 0.04*1e6;
    R = [0.15 0.50 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 200;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -6000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 14
    ntot = 0.04*1e6;
    R = [0.15 0.50 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 200;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -5000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 15
    ntot = 0.04*1e6;
    R = [0.25 0.40 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 200;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -5000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 16
    ntot = 0.04*1e6;
    R = [0.25 0.40 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -4000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 17
    ntot = 0.04*1e6;
    R = [0.35 0.30 0.35]; 
    n = ntot.*R;        
    
    T1 = 30;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -4000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 18 % new 2021-08-30
    ntot = 0.04*1e6;
    R = [0.65 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -7000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd3];
    T = [T1 T3];
  case 19 % new 2021-08-30
    ntot = 0.04*1e6;
    R = [0.4 0.25 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -2000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
  case 20 % new 2021-08-30
    ntot = 0.04*1e6;
    R = [0.4 0.25 0.35]; 
    n = ntot.*R;        
    
    T1 = 50;
    T2 = 100;
    T3 = 1500;
    vd1 = -11000*1e3;
    vd2 = -2000*1e3;
    vd3 = 6000*1e3;  
    
    n = n;
    vd = [vd1 vd2 vd3];
    T = [T1 T2 T3];
end


units =irf_units;
vt = sqrt(2*units.e*T./units.me); % m/s

nPop = numel(n);
if 0
f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
else
f0_str = ['f0 = @(v,n,vd,vt) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
end
eval(f0_str)

params.n = n;
params.vd = vd;
params.T = T;
params.vt = vt;
params.R = R;