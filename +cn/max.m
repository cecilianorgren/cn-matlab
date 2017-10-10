function out = max(A)
maxInd = find(abs(A)==max(abs(A)));
out = A(maxInd);
out = max(out); % choose the positive if they have exact same value