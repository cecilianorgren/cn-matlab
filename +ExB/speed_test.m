n = 1e5;

tic
x = randn(n,1);
for k = 1:n
    ExB.speed_test_function(x(k),x(k),x(k));
end
toc % Results: 3.0
tic
for k = 1:n
    ExB.speed_test_function(0.5,0.5,0.5);
end
toc % Results: 2.7
tic
for k = 1:n
    ExB.speed_test_function;
end
toc % Results: 2.4  

