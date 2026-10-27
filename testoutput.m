function out = testoutput()

PD = TSeries(EpochTT(1:10),(1:10)');
out = {PD,PD*2};
