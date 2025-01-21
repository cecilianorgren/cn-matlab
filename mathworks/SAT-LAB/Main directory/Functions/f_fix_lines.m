function [phi_corr,lamda_corr] = f_fix_lines(phi,lamda,ang_threshold)
%%
% F_FIX_LINES corrects the west-east line bug that appears when a line
% passes from eastern hemisphere to western hemisphere or the opposite.
% Longitude coordinated should be from 0 to 360.
%
% HOW: [phi_corr,lamda_corr] = f_fix_lines(phi,lamda,ang_threshold)
%
% Input: phi            [1 x n] latitude in degrees.
%                    or [n x 1]
%
%        lamda          [1 x n] longitude in degrees.
%                    or [n x 1]
%
%        ang_threshold  [1 x 1] angle difference threshold in degrees to
%                               apply the correction.
%
% Output: phi_corr      [1 x n] corrected latitude in degrees.
%                    or [n x 1]
%
%         lamda_corr    [1 x n] longitude in degrees.
%                    or [n x 1]
%
% Dimitrios Piretzidis, Department of Geomatics Engineering, UofC
% 23/01/2017
%
% uses m-files: none

%% Revision history

%% Remarks

%% Start the algorithm

%Input check
if nargin ~= 3
    error('Wrong number of input arguments')
end

if max(size(phi) ~= size(lamda) ) == 1
    error('Input coordinates do not have the same dimensions')
end

if max(size(ang_threshold)) ~= 1
    error('Select only one angle threshold')
end

phi                                   = phi(:);
lamda                                 = lamda(:);

diff_lamda                            = diff(lamda);
[rr,~]                                = find(abs(diff_lamda) > ang_threshold);
rr                                    = rr + 1;

phi_corr                              = nan(size(phi,1) + 3*size(rr,1),1);
lamda_corr                            = nan(size(phi,1) + 3*size(rr,1),1);
su                                    = 0;

for kk = 1:size(phi,1)
    
    if max(kk == rr) == 1
        
        if diff_lamda(kk-1,1) < ang_threshold
            
            phi_corr(kk + su,1)       = phi(kk,1);
            
            if max(lamda_corr) > 180
                
                lamda_corr(kk + su,1) = 359.9999;
                
            end
            
        elseif diff_lamda(kk-1,1) > ang_threshold
            
            phi_corr(kk + su,1)       = phi(kk,1);
            
            if max(lamda_corr) > 180
                
                lamda_corr(kk + su,1) = 0;
                
            end
            
        end
        
        su                            = su + 1;
        
        su                            = su + 1;
        
        if diff_lamda(kk-1,1) < ang_threshold
            
            phi_corr(kk + su,1)       = phi(kk,1);
            
            if max(lamda_corr) > 180
                
                lamda_corr(kk + su,1) = 0;
                
            end
            
        elseif diff_lamda(kk-1,1) > ang_threshold
            
            phi_corr(kk + su,1)       = phi(kk,1);
            
            if max(lamda_corr) > 180
                
                lamda_corr(kk + su,1) = 359.9999;
                
            end
            
        end
        
        su                            = su + 1;
        
    end
    
    phi_corr(kk + su,1)               = phi(kk,1);
    lamda_corr(kk + su,1)             = lamda(kk,1);
    
end

if size(phi,1) == 1
    
    phi_corr                          = phi_corr';
    lamda_corr                        = lamda_corr';
    
end

end
