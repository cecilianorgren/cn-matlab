function out = mEV(pot,range,equity,mV)

mEV = (1-range)*pot+range*(equity*pot+(1-equity)*mV);