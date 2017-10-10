function out = abcd(a,b,c,d)

ab=a+b;
out =c+ab;
if ~isempty(d)
    out=out+d;
end