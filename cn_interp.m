function val = cn_interp(t,ts,flag)


k=diff(ts(:,2))/diff(ts(:,1));
dt=ts(end,1)-ts(1,1);
if flag==1
    val=ts(1,2)+k*(t-ts(1,1));
else
    val=[t ts(1,2)+k*(t-ts(1,1))];
end
