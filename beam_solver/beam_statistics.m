% beam_statistics
% script running some beam solver and extracts vph vT etc

k = 0.01:0.05:2;
n = [1 0.5 0.5]; 
T = [1000 1000 600];
vd = [0 -1.7*cn_eV2v(1000,'eV') 1.7*cn_eV2v(1000,'eV')];
TolFun = 1e-8;
doPlot = 1;
[wr_out,wi_out,wr_allout,wi_allout,wrin_out,wiin_out,param_out] = beam_solver_bistream(n,T,vd,k,doPlot,TolFun);

