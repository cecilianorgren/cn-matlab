% input data
Te1 = 1600;
Te2 = [100 40];
Ti = [800 1600];
S = 0.1:0.05:1.5;
R = [0.05:0.95 0.99];
k = 0.01:0.02:2.2; 
doPlot = 0;
%%
if 0    
    % input data
    Te1 = 1600;
    Te2 = [60];
    Ti = [1600];
    S = 0.2:0.1:1.2;
    R = [0.19:0.1:0.99];
    k = 0.01:0.05:0.8; 
    doPlot = 0;
end
if 0
    % input data
    Te1 = 1600;
    Te2 = [60];
    Ti = [1600];
    S = 0.3:0.1:1;
    R = [0.3 0.7];
    k = 0.01:0.1:1.5; 
    doPlot = 0;
end
if 1
    % input data
    Te1 = 1600;
    Te2 = [60];
    Ti = [1600];
    S = 0.3:0.1:1.4;
    R = [0.2 0.3];
    k = 0.01:0.1:2; 
    doPlot = 0;
end

if 1
    % input data
    Te1 = 300;
    Te2 = Te1/25;
    Ti = Te1;
    S = 0.1:0.05:1.5;
    R = [0.05:0.05:0.95 0.99];
    k = 0.01:0.02:2.2; 
    doPlot = 0;
    n = 1;
end

% run the solver
%n = 0.06;
TolFun = 1e-8;
%TolFun = 1e-12;
tic
[wr,wi,wrall,wiall,wrin,wiin,wrinall,wiinall,param] = beam_solver(Te1,Te2/Te1,Ti/Te1,R,S,k,n,doPlot,TolFun);
toc
[kmax,wrmax,wimax,vphmax,vph] = beam_solver_max(k,wr,wi,0);
disp('Done!')
clear doPlot;
nTe1 = numel(Te1);
nTe2 = numel(Te2);
nTi = numel(Ti);
%% save
savePath = '/home/cecilia/Research/BeamSolver/';
saveName = [datestr(now,'yyyy-mm-ddTHHMMSS') '_beam_solver_' num2str(numel(k)) 'kx' num2str(numel(S)) 'Sx' num2str(numel(R)) 'Rx' num2str(numel(Te2)) 'Te2x' num2str(numel(Ti)) 'Ti'];
save([savePath saveName])
disp(['Saved ' savePath saveName '.mat'])
