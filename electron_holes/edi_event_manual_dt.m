function manual = edi_event_manual_dt()
% Speeds are in locally parallel coordinates dfined by perp1, perp2, par.
% tref is approximately the center of the ESW as seen by mms1, so that
% t1 = t_ref + dt(1), etc
% tminus and tplus gives approximate regions where I could perhaps read off
% the potential. It's not altogether clear where these ponits should be
% tp1, tp2, approximate times of first and second peak.
ieh = 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.533947021Z';
manual(ieh).tminus = '2017-07-06T13:54:05.532483886Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.535222412Z'; 
manual(ieh).tp1 = [-0.00050,-0.00033,-0.00023,-0.00031];
manual(ieh).tp2 = [ 0.00045, 0.00031, 0.00025, 0.00035];
manual(ieh).gseBav = [-22.7878    1.6498   -1.2820];
manual(ieh).perp1 = [0.0720    0.9974    0.0040];
manual(ieh).perp2 = [0.0562         0   -0.9984];
manual(ieh).par   = [-0.9958    0.0721   -0.0560];  
manual(ieh).dt  = [0.0000,-0.00115,-0.00087,-0.00132]; 
manual(ieh).v  = 8.78e+03*[1.00,-0.04, 0.05];

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.538696044Z'; 
manual(ieh).tminus = '2017-07-06T13:54:05.537017822Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.540224609Z'; 
manual(ieh).tp1 = [-0.00076,-0.00051,-0.00054,-0.00058];
manual(ieh).tp2 = [ 0.00043, 0.00048, 0.00030, 0.00037];
manual(ieh).gseBav = [-22.7882    1.6481   -1.2911];
manual(ieh).perp1 = [0.0719    0.9974    0.0041];
manual(ieh).perp2 = [0.0566   -0.0000   -0.9984];
manual(ieh).par   = [-0.9958    0.0720   -0.0564];   
manual(ieh).dt  = [0.0000   -0.00075   -0.0011   -0.0015]; 
manual(ieh).v = 8.02e+03 * [0.94 -0.35 -0.03]; % 

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.547171386Z';
manual(ieh).tminus = '2017-07-06T13:54:05.545546875Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.548611816Z';
manual(ieh).tp1 = [-0.00049,-0.00049,-0.00051,-0.00051];
manual(ieh).tp2 = [ 0.00050, 0.00059, 0.00063, 0.00059];
manual(ieh).gseBav = [-22.7855    1.6436   -1.3151];
manual(ieh).perp1 = [0.0717    0.9974    0.0041];
manual(ieh).perp2 = [0.0576    0.0000   -0.9983];
manual(ieh).par   = [-0.9958    0.0718   -0.0575];   
manual(ieh).dt  = [0.0000   -0.00162   -0.00108   -0.00125]; 
manual(ieh).v = 7.62e+03 * [0.98  0.20 -0.05]; % this one is ok, clearly something wrong with EH_properties7.62e+03 * [0.98  0.20 -0.05]

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.549964843Z';
manual(ieh).tminus = '2017-07-06T13:54:05.548452392Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.551060302Z';
manual(ieh).tp1 = [-0.00034,-0.00032,-0.00034,-0.00040];
manual(ieh).tp2 = [ 0.00034, 0.00025, 0.00041, 0.00034];
manual(ieh).gseBav = [-22.7882    1.6481   -1.2911];
manual(ieh).perp1 = [0.0718    0.9974    0.0042];
manual(ieh).perp2 = [0.0579         0   -0.9983];
manual(ieh).par   = [-0.9957    0.0720   -0.0578]; 
manual(ieh).gseBav = [ -22.7852    1.6469   -1.3216]; 
manual(ieh).dt  = [0.0000   -0.0013   -0.00093   -0.0013];
manual(ieh).v = 8.44e+03 * [1.00  0.05  0.02]; % 

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.556753417Z';
manual(ieh).tminus = '2017-07-06T13:54:05.555059570Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.558725585Z';
manual(ieh).tp1 = [-0.00042,-0.00040,-0.00040,-0.00026];
manual(ieh).tp2 = [ 0.00067, 0.00048, 0.00046, 0.00036];
manual(ieh).gseBav = [-22.7891    1.6627   -1.3376];
manual(ieh).perp1 = [0.0725    0.9974    0.0043];
manual(ieh).perp2 = [0.0586   -0.0000   -0.9983];
manual(ieh).par   = [-0.9956    0.0726   -0.0584]; 
manual(ieh).gseBav = [ -22.7891    1.6627   -1.3376]; 
manual(ieh).dt  = [0.0000   -0.00105   -0.00115   -0.0012]; 
manual(ieh).v = 8.74e+03 * [0.98 -0.09 -0.15]; % 

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.562840820Z';
manual(ieh).tminus = '2017-07-06T13:54:05.561382324Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.564882080Z';
manual(ieh).tp1 = [-0.00054,-0.00047,-0.00050,-0.00063];
manual(ieh).tp2 = [ 0.00034, 0.00046, 0.00032, 0.00046];
manual(ieh).gseBav = [-22.7938    1.6766   -1.3445];
manual(ieh).perp1 = [0.0731    0.9973    0.0043];
manual(ieh).perp2 = [0.0589         0   -0.9983];
manual(ieh).par   = [-0.9956    0.0732   -0.0587]; 
manual(ieh).dt  = [0.0000,-0.00135,-0.00098,-0.00145];  
manual(ieh).v  = 7.81e+03*[1.00  0.00  0.05];

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.567746826Z';
manual(ieh).tminus = '2017-07-06T13:54:05.566484863Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.569294921Z';
manual(ieh).tp1 = [-0.00037,-0.00048,-0.00024,-0.00042];
manual(ieh).tp2 = [ 0.00043, 0.00043, 0.00041, 0.00072];
manual(ieh).gseBav = [-22.8033    1.6873   -1.3534];
manual(ieh).perp1 = [0.0735    0.9973    0.0044];
manual(ieh).perp2 = [0.0592   -0.0000   -0.9982];
manual(ieh).par   = [-0.9955    0.0737   -0.0591]; 
manual(ieh).dt  = [0.0000,-0.00063,-0.00080,-0.0014];  
manual(ieh).v  = 9.12e+03 * [0.92 -0.37  0.09];

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.572459472Z';
manual(ieh).tminus = '2017-07-06T13:54:05.571583496Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.573349365Z';
manual(ieh).tp1 = [-0.00025,-0.00061,-0.00045,-0.00025];
manual(ieh).tp2 = [ 0.00025, 0.00071, 0.00030, 0.00023];
manual(ieh).gseBav = [-22.8114    1.6932   -1.3658];
manual(ieh).perp1 = [0.0738    0.9973    0.0044];
manual(ieh).perp2 = [0.0598         0   -0.9982];
manual(ieh).par   = [-0.9955    0.0739   -0.0596]; 
manual(ieh).dt  = [0.0000,-0.00075,-0.00120,-0.0013];  
manual(ieh).v  = 8.51e+03*[0.94,-0.30,-0.16];

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.579765136Z';
manual(ieh).tminus = '2017-07-06T13:54:05.578395263Z'; 
manual(ieh).tplus = '2017-07-06T13:54:05.581858886Z';
manual(ieh).tp1 = [-0.00056,-0.00024,-0.00051,-0.00047];
manual(ieh).tp2 = [ 0.00064, 0.00072, 0.00065, 0.00055];
manual(ieh).gseBav = [-22.8235    1.7000   -1.3811];
manual(ieh).perp1 = [0.0740    0.9972    0.0045];
manual(ieh).perp2 = [0.0604         0   -0.9982];
manual(ieh).par   = [-0.9954    0.0741   -0.0602]; 
manual(ieh).dt  = [0.0000,-0.00097,-0.00086,-0.0012];  
manual(ieh).v  = 9.66e+03*[1.00,-0.10, 0.01];

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.585509765Z';
manual(ieh).tminus = '2017-07-06T13:54:05.583967773Z';
manual(ieh).tplus = '2017-07-06T13:54:05.586859374Z';
manual(ieh).tp1 = [-0.00037,-0.00034,-0.00036,-0.00039];
manual(ieh).tp2 = [ 0.00037, 0.00034, 0.00039, 0.00050];
manual(ieh).gseBav = [-22.8265    1.7105   -1.3870];
manual(ieh).perp1 = [0.0745    0.9972    0.0045];
manual(ieh).perp2 = [0.0606   -0.0000   -0.9982];
manual(ieh).par   = [-0.9954    0.0746   -0.0605]; 
manual(ieh).dt = [0.0000,-0.00085,-0.00115,-0.00145]; 
manual(ieh).v = 8.06e+03*[0.96,-0.28,-0.07];

ieh = ieh + 1;
manual(ieh).t_ref ='2017-07-06T13:54:05.604416748Z';
manual(ieh).tminus = '2017-07-06T13:54:05.603427490Z';
manual(ieh).tplus = '2017-07-06T13:54:05.605484374Z';
manual(ieh).tp1 = [-0.00023,-0.00029,-0.00016,-0.00029];
manual(ieh).tp2 = [ 0.00025, 0.00037, 0.00024, 0.00025];
manual(ieh).gseBav = [ -22.8326    1.7458   -1.4273];
manual(ieh).perp1 = [0.0759    0.9971    0.0047];
manual(ieh).perp2 = [0.0624   -0.0000   -0.9981];
manual(ieh).par   = [-0.9952    0.0761   -0.0622]; 
manual(ieh).dt = [0.0000,-0.00095,-0.00075,-0.0011]; 
manual(ieh).v = 1.05e+04*[1.00,-0.05, 0.04];

%   ieh = ieh + 1;
%   manual(ieh).t_ref = '2017-07-06T13:54:05.615465576Z';
%   manual(ieh).gseBav = [-22.8243    1.7780   -1.4472];
%   manual(ieh).perp1 = [0.0774    0.9970    0.0049];
%   manual(ieh).perp2 = [0.0633    0.0000   -0.9980];
%   manual(ieh).par   = [-0.9950    0.0775   -0.0631]; 
%   manual(ieh).dt = [nan nan nan nan]; 
%   manual(ieh).v = nan * [nan nan nan]; % unclear timing here, only ok between 3 and 4

ieh = ieh + 1;
manual(ieh).t_ref = '2017-07-06T13:54:05.672972656Z';
manual(ieh).tminus = '2017-07-06T13:54:05.671991699Z';
manual(ieh).tplus = '2017-07-06T13:54:05.674228759Z';
manual(ieh).tp1 = [-0.00031,-0.00020,-0.00024,-0.00038];
manual(ieh).tp2 = [ 0.00026, 0.00022, 0.00036, 0.00033];
manual(ieh).gseBav = [-22.8477    1.8595   -1.5054];
manual(ieh).perp1 = [0.0808    0.9967    0.0053];
manual(ieh).perp2 = [0.0657         0   -0.9978];
manual(ieh).par   = [-0.9946    0.0809   -0.0655]; 
manual(ieh).dt = [0.0000,-0.00083,-0.00095,-0.00123]; 
manual(ieh).v = 9.48e+03*[0.98,-0.21,-0.05];

  % Calculate and add a few extra paramters
for ieh = 1:numel(manual)  
  tref = EpochTT(manual(ieh).t_ref);
  dt = manual(ieh).dt;
  tp1 = tref + manual(ieh).tp1;
  tp2 = tref + manual(ieh).tp2;
  tpp = tp2-tp1;
  lmn = [manual(ieh).perp1; manual(ieh).perp2; manual(ieh).par]; 
  gsev = manual(ieh).v; % velocity is given in gse coordinate system
  pppv = gsev*lmn';
  pppvnorm = pppv/norm(pppv); % unit vector of velocity
  vpar = pppv(3);    
  pppBav = manual(ieh).gseBav*lmn';
  
  manual(ieh).pppB = pppBav;
  manual(ieh).vpar = vpar;
  manual(ieh).tpp = tpp;
  manual(ieh).lpp = tpp*vpar;    
  manual(ieh).pitchangle = acosd(pppvnorm(3));
end