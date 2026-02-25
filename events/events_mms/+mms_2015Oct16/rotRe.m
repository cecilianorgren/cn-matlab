function out = gradP(R1,R2,R3,R4,Re1,Re2,Re3,Re4)
% R are TSeries
% P are 3x3 cell array of TSeries E+VexB

ts = Re1.time;
c_eval('R?=R?.resample(ts);')
c_eval('Re?=Re?.resample(ts);',2:4)
rotRe = c_4_grad('R?','Re?','curl');
out = rotRe;