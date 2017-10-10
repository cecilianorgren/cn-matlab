pa1=15:30:180;
pa2=7.5:15:180;
a=[0 90 180];
na=length(a);
for k=1:na
    a(k)
    diff=abs(pa1-a(k))
    ind=find(diff==min(diff))
end