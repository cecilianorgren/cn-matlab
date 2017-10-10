% define k 
k = 0:0.01:2;
wr = [0];
wi = [0];


% zero values



af = @(temp) beam_solver_disprel(temp,k(ik),fv);
