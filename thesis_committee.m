rng default
A = randn(100000,1);
iqrA = iqr(A);
hist(A,30)