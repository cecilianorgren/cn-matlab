% XES.EVENTS    Input parameters to xes1 runs.
%   

% Two-stream script 
% Variable parameters

% Physical variables:
event = 'R0.30 S0.50';
no = 0;
doInput = 1;
switch event
    case '1' % Daniels fast ESWs, original
        B = 25;                               
        n1 = 0.01; n2 = 0.003; n = n1 + n2;
        R = n2/n;
        Te1 = 4000; Te2  = 250; Ti = 2000;
        vb_in_eV = 1200;
        S = sqrt(vb_in_eV/Te1);
    case '1b' % Daniels fast ESWs, lower beam temp
        B = 25;                               
        n1 = 0.01; n2 = 0.003; n = n1 + n2;
        R = n2/n;
        Te1 = 4000; Te2  = 100; Ti = 2000;
        vb_in_eV = 1200;
        S = sqrt(vb_in_eV/Te1);
    case '1c' % Daniels fast ESWs, lower beam temp and higher beam speed
        B = 25;                               
        n1 = 0.01; n2 = 0.003; n = n1 + n2;
        R = n2/n;
        Te1 = 4000; Te2  = 100; Ti = 2000;
        vb_in_eV = 2000;
        S = sqrt(vb_in_eV/Te1);
    case '1d' % Daniels fast ESWs, lower beam temp and higher beam speed, warmer background
        B = 25;                               
        n1 = 0.01; n2 = 0.003; n = n1 + n2;
        R = n2/n;
        Te1 = 8000; Te2  = 100; Ti = 2000;
        vb_in_eV = 2000;
        S = sqrt(vb_in_eV/Te1);
    case '2' % Daniels slow ESWs, original
        B = 25;
        vb_in_eV = 200;               
        n1 = 0.12; n2 = 0.03; n = n1 + n2; R = n2/n;
        Te1 = 3500; Te2  = 120; Ti = 5000;
        S = sqrt(vb_in_eV/Te1);
        S = cn_eV2v(vb_in_eV,'eV')/cn_eV2v(Te1,'eV');
    case '2b' % Daniels slow ESWs
        B = 25;
        vb_in_eV = 300;               
        n1 = 0.12; n2 = 0.03; n = n1 + n2; R = n2/n;
        Te1 = 3500; Te2  = 60; Ti = 5000;
        S = sqrt(vb_in_eV/Te1);
        S = cn_eV2v(vb_in_eV,'eV')/cn_eV2v(Te1,'eV');
    case '3' 
        S=1.6; % beam drift evlocity in units of background thermal velocity
        B = 25;
        n = 0.06;
        R = 0.99;
        Te1 = 1600;
        Te2 = 60;
        Ti = 1600;
        n1= n*(1-R);
        n2= n*R;        
    case '4' % Daniels slow ESWs, lower beam temperature
        B = 25;
        vb_in_eV = 200;        
        n = 0.15;        
        n1 = 0.12;
        n2 = 0.03;
        R = n2/n;
        Te1 = 3500; 
        Te2  = 20; 
        Ti = 5000;
        S = sqrt(vb_in_eV/Te1);        
    case '5' % Daniels fast ESWs, lower beam temp
        B = 25;
        vb_in_eV = 1200;
        S = 0.2;
        n = 0.1;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 100; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);   
    case '6' % Daniels fast ESWs, even lower beam temp
        B = 25;
        vb_in_eV = 1200;
        S = 0.2;
        n = 0.1;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);   
    case '7' % Daniels fast ESWs, even lower beam temp
        B = 25;
        vb_in_eV = 2000;
        S = 0.2;
        n = 0.1;        
        n1 = 0.01;
        n2 = 0.003;
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);  
    case '8' % Daniels fast ESWs, even lower beam temp and lower density
        B = 25;
        vb_in_eV = 3000;
        S = 0.2;               
        n1 = 0.01;
        n2 = 0.001;
        n = n1+n2; 
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1); 
    case '9' % Daniels fast ESWs, even lower beam temp and lower density
        B = 25;
        vb_in_eV = 3000;
        S = 0.2;               
        n1 = 0.009;
        n2 = 0.001;
        n = n1+n2; 
        R = n2/n;
        Te1 = 4000; 
        Te2  = 50; 
        Ti = 2000;
        S = sqrt(vb_in_eV/Te1);  
    case 'mp1' % Daniels slow and fast ESWs, magnetopause
        B = 25;                               
        n1 = 1.8; n2 = 0.4; n = n1 + n2;
        R = n2/n;
        Te1 = 110; Te2  = 12; Ti = 220;
        vb_in_eV = 90;
        S = sqrt(vb_in_eV/Te1);
    case 'mp2' % Daniels slow and fast ESWs, magnetopause, higher Tbg       
        B = 25;                               
        n1 = 1.8; n2 = 0.4; n = n1 + n2;
        R = n2/n;
        Te1 = 310; Te2  = 12; Ti = 2*Te1;
        vb_in_eV = 90;
        S = sqrt(vb_in_eV/Te1);    
    case 'mp3' % Daniels new MP event, slow ESWs, include in paper
        B = 25;                               
        n1 = 0.2+0.1; n2 = 0.1; n = n1 + n2; % observed is n2 = 0.3
        R = n2/n;
        Te1 = 350; Te2  = 14; Ti = 2*Te1; % observed Te2 is 80
        vb_in_eV = 90;
        S = sqrt(vb_in_eV/Te1); 
    case 'mp3b' % Daniels new MP event, slow ESWs, include in paper
        B = 25;                               
        n1 = 0.2+0.1; n2 = 0.2; n = n1 + n2; % observed is n2 = 0.3
        R = n2/n;
        Te1 = 350; Te2  = 14; Ti = 2*Te1; % observed Te2 is 80
        vb_in_eV = 90;
        S = sqrt(vb_in_eV/Te1); 
    case 'mp3c' % Daniels new MP event, slow ESWs, include in paper
        B = 25;                               
        n1 = 0.2+0.1; n2 = 0.15; n = n1 + n2; % observed is n2 = 0.3
        R = n2/n;
        Te1 = 350; Te2  = 14; Ti = 2*Te1; % observed Te2 is 80
        vb_in_eV = 90;
        S = sqrt(vb_in_eV/Te1); 
    case 'R0.20 S0.50' % Daniels new MP event, slow ESWs, include in paper
        n=1;
        R=0.2;
        S=0.5;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1; % observed Te2 is 80                
    case 'R0.25 S0.50' % 
        n=1;
        R=0.25;
        S=0.5;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    case 'R0.30 S0.50' % 
        n=1;
        R=0.3;
        S=0.5;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    case 'R0.25 S0.40' % 
        n=1;
        R=0.25;
        S=0.4;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    case 'R0.25 S0.45' % 
        n=1;
        R=0.25;
        S=0.45;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    case 'R0.20 S0.40' % 
        n=1;
        R=0.20;
        S=0.4;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    case 'R0.20 S0.35' % 
        n=1;
        R=0.20;
        S=0.35;
        n1 = n*(1-R);
        n2 = n*R;
        Te1 = 300; Te2  = 12; Ti = Te1;
    otherwise
        disp('Did not recognize case name.')
        doInput = 0;
end
R
if doInput
%try
    disp(['--- ' event ' ---'])
    xes.input([n1 n2 n1+n2],[Te1 Te2 Ti],[0 S*cn_eV2v(Te1,'eV') 0],{'e','e','i'})
%end
end