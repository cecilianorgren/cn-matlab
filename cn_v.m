function [v vminus vplus] = cn_v(t0,terror,vhat,Pos3,Pos4,t)
c_eval('gsePos?=cn_toepoch(t,Pos?);',3:4);
r34=gsePos4(2:4)-gsePos3(2:4); % Relative position
v=cn_scalar(r34,vhat)/t0;
vminus=cn_scalar(r34,vhat)/(t0+terror);
vplus=cn_scalar(r34,vhat)/(t0-terror);