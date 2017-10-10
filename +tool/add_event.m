% Load TimeTable
load TT
% includes movingTT and singleTT
% movingTT includes longer time intervals
% singleTT includes shorter time intervals

switch toSave
    case 'single'
        disp('Adding event to singleTT.')
        clear singleDescription
        singleDescription.SC = ['C' num2str(sc)];
        singleDescription.tint = tint;
        singleDescription.vecC = corr_dir;
        singleDescription.maxC = max(corr_dir(:,1));
        singleDescription.k = direction;
        singleDescription.B0 = z(1,:)*B0;        
        singleDescription.v = velocity; % Velocity that minimizes area between phi_E and phi_B
        singleDescription.phiE = phi_E(:,[1 1+i_v]); % Potential at 'best' v
        singleDescription.phiB = phi_B;
        singleDescription.flh = flh_loc;
        singleDescription.ffilt = f_highpass;
        singleDescription.mva_l = mva_l;
        singleDescription.mva_v = mva_v;
        singleDescription.angles = angles;
        singleDescription.n = n_loc;
        
        singleComment = comment;
        singleTT = add(singleTT,tint,singleDescription,singleComment);
    case 'moving'
        movingDescription.tint = tint;
        %movingDescription.
        %movingDescription.
        %movingDescription.
        %movingDescription.
        %movingDescription.
        movingTT = add(movingTT,tint,description,comment);
end

% Save TimeTable again
save('TT','movingTT','singleTT')