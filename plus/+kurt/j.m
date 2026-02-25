function [j jrot] = j(r1,r2,r3,r4,b1,b2,b3,b4,angledeg) 
% Estimate the current from a four position spacecraft.
% Due to coordinate system construction, it is necessary to first rotate
% all the vectors (to avoid division by zero), and then rotate them back.
% Let's se if it works.

%angle = pi/4; % 45 degrees
angle = angledeg*2*pi/360;

c_eval('r?(:,2) = cos(angle)*r?(:,2) + sin(angle)*r?(:,3);',3:4)
c_eval('r?(:,3) = cos(angle)*r?(:,3) - sin(angle)*r?(:,2);',3:4)

c_eval('b?(:,2) = cos(angle)*b?(:,2) + sin(angle)*b?(:,3);',3:4)
c_eval('b?(:,3) = cos(angle)*b?(:,3) - sin(angle)*b?(:,2);',3:4)

[jrot,divB,B,jxB,divTshear,divPb] = c_4_j(r1,r2,r3,r4,b1,b2,b3,b4);

j = jrot;
c_eval('j(:,2) = cos(-angle)*jrot(:,2) + sin(-angle)*jrot(:,3);',3:4)
c_eval('j(:,3) = cos(-angle)*jrot(:,3) - sin(-angle)*jrot(:,2);',3:4)





