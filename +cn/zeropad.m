function out = zeropad(A,tint)
%CN.ZEROPAD    Zeropads.
%   CN.ZEROPAD(A,tint) zeropads the timeseries A out to the interval tint
%   with the same sampling frequency as A.

f = cn.f(A);
dt = 1/f;

leftTs = tocolumn(tint(1):dt:(A(1,1)-dt));
rightTs = tocolumn((A(end,1)+dt):dt:tint(2));

leftZeros = [leftTs leftTs*0];
rightZeros = [rightTs rightTs*0];

out = [leftZeros; A; rightZeros];

    