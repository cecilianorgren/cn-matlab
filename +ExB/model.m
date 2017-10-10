switch modelnr
case 1 % 2007-08-31        
    name = '2007-08-31'; lr   = 15; lz   = 10;
    phi0 = 300; % changed this from 200 to 300  2014-03-24
    B0   = 25; Tpar = 1600; Tper = 1600; veh  = 500; n = 0.07; % cc
case 2; name = 'Tao2011'; 
    lr   = 35; lz   = 35; phi0 = 3000;
    B0   = 50; Tpar = 8000; Tper = 5000; veh  = 100000; n = 0.02; % cc
case 3; name = 'Basic'; 
    lr   = 20; lz   = 20; phi0 = 0.01; B0   = 60; Tpar = 1000;
    Tper = 1000; veh  = 0; n = 0.05;
    % Tper = 1000; B0 = 10 gives rho_e ~ 10.        
case 4; name = 'ExagTao2011';
    lr = 35; lz = 35; phi0 = 10000; B0 = 50; Tpar = 8000; Tper = 5000; veh = 0; n = 0.02;
case 5; name = 'WrongTao';
    lr = 3.5; lz = 3.5; phi0 = 300; B0 = 50; Tpar = 8000; Tper = 5000; veh = 10000; n = 0.02;
case 6; name = 'WrongTao_half_speed';
    lr = 18; lz = 18; phi0 = 1500; B0 = 50; Tpar = 8000; Tper = 5000; veh = 50000;  n = 0.02;
case 7; name = '2007-08-31_tweak';
    lr = 15; lz = 15; phi0 = 1000; B0 = 25; Tpar = 1600; Tper = 1600; 
    veh  = 500; n = 0.07; % cc            
case 8; name = '2007-08-31_biggersize';
    lr = 30; lz = 20; phi0 = 600; B0 = 25; Tpar = 1600; Tper = 1600; 
    veh  = 500; n = 0.07/2; % cc                    
case 9; name = '2007-08-31_coolerperp';
    lr = 15; lz = 10; phi0 = 300; B0 = 25; 
    Tpar = 1600; Tper = 1600/2; veh  = 500; n = 0.07; % cc 
case 10; name = '2007-08-31_hotterperp';
    lr = 15; lz = 10; phi0 = 300; B0 = 25;
    Tpar = 1600; Tper = 1600*4; veh = 500; n = 0.07; % cc   
case 11; name = '2007-08-31_cooler';
    lr = 15; lz = 10; phi0 = 300; B0 = 25; 
    Tpar = 1600/2; Tper = 1600/2; veh  = 500; n = 0.07; % cc 
case 12; name = '2007-08-31_hotter';
    lr = 15; lz = 10; phi0 = 300; B0   = 25;
    Tpar = 1600*4; Tper = 1600*4; veh  = 500; n = 0.07; % cc       
case 13; name = '2007-08-31_hotterpar';
    lr   = 15; lz   = 10; phi0 = 300; B0   = 25;
    Tpar = 1600*8; Tper = 1600; veh  = 500; n = 0.07; % cc       
case 14; name = '2007-08-31_hotterperp';
    lr   = 15; lz   = 10; phi0 = 300; B0   = 25;
    Tpar = 1600; Tper = 1600*2; veh  = 500; n = 0.07; % cc     
case 15; name = '2007-08-31_coolerpar';
    lr   = 15; lz   = 10; phi0 = 300; B0   = 25;
    Tpar = 1600/2; Tper = 1600; veh  = 500; n = 0.07; % cc       
case 16 % 2007-08-31        
    name = '2007-08-31-0V';
    lr   = 15;
    lz   = 10;
    phi0 = 0;
    B0   = 25;
    Tpar = 1600;
    Tper = 1600;        
    veh  = 500;
    n = 0.07; % cc
case 17 % Tao2011
    name = 'Tao2011-0V';
    lr   = 35;
    lz   = 35;
    phi0 = 0;
    B0   = 50;
    Tpar = 8000;
    Tper = 5000;        
    veh  = 100000;
    n = 0.02; % cc
case 18; name = '2007-08-31_hotterpar2';
    lr   = 15; lz   = 10; phi0 = 300; B0   = 25;
    Tpar = 1600*16; Tper = 1600; veh  = 500; n = 0.07; % cc 
case 19; name = '2007-08-31-wider-less-dense';
    lr   = 25; lz   = 10; phi0 = 300; B0 = 25;
    Tpar = 1600; Tper = 1600; veh  = 500; n = 0.035; % cc  
case 20; name = '2007-08-31-200V';
    lr   = 15; lz   = 10; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 1600; veh  = 500; n = 0.035; % cc  
case 21; name = '2007-08-31-200V_lr0_lz5';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 1600; veh  = 500; n = 0.035; % cc  
case 22; name = '2007-08-31-200V_lr0_lz5_16eV';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 16; veh  = 500; n = 0.035; % cc  
case 23; name = '2007-08-31-200V_lr0_lz5_160eV';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 160; veh  = 500; n = 0.035; % cc  
case 24; name = '2007-08-31-200V_lr0_lz5_500eV';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 500; veh  = 500; n = 0.035; % cc  
case 25; name = '2007-08-31-200V_lr0_lz5_1300eV';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 1300; veh  = 500; n = 0.035; % cc  
case 26; name = '2007-08-31-200V_lr0_lz5_2500eV';
    lr   = 9; lz   = 5; phi0 = 200; B0 = 25;
    Tpar = 1600; Tper = 2500; veh  = 500; n = 0.035; % cc          
case 27; name = '2007-08-31-500V_lr12_lz5_2000eV';
    lr   = 12; lz   = 5; phi0 = 500; B0 = 25;
    Tpar = 2000; Tper = 2000; veh  = 500; n = 0.035; % cc          
case 28; name = '2007-08-31-500V_1600eV_9_5';
    lr   = 12; lz   = 5; phi0 = 500; B0 = 25;
    Tpar = 1600; Tper = 1600; veh  = 500; n = 0.035; % cc  
end     

if exist('printModel','var') 
    if printModel        
        mu0 = 1.2566e-6; e = 1.6022e-19;
        tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT
        disp([num2str(modelnr) '. ' name ': tdB=' num2str(tdB) ', lr=' num2str(lr) ,...
            ', lz=' num2str(lz) ', phi0=' num2str(phi0) ', B0=' num2str(B0),...
            ', Tpar=' num2str(Tpar) ', Tper=' num2str(Tper) ', veh=' num2str(veh) ', n=' num2str(n)])
    end
end