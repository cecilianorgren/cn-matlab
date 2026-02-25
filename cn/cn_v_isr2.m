function [vx vy] = cn_v_isr2(t,v,vd,way)
x=c_coord_trans('DSI','GSE',[toepoch(t(1:6)) [1 0 0]],'CL_ID',3);
y=c_coord_trans('DSI','GSE',[toepoch(t(1:6)) [0 1 0]],'CL_ID',3);

cosalfax=cn_scalar(x(2:4),vd); % cos alfa
cosalfay=cn_scalar(y(2:4),vd);

switch way
    case 1
        vx=v/cosalfax;
        vy=v/cosalfay;
    case 2
        vx=v*cosalfax;
        vy=v*cosalfay;
end