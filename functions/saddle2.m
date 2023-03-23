function indices = saddle2(M)
saddle_m =  zeros(size(M));
row_maxima = max(M, [], 2);
col_minima = min(M, [], 1);
indices = [];
for i = 1:size(M, 1)
    for j = 1:size(M, 2)
        if M(i, j) == row_maxima(i) && M(i, j) == col_minima(j)
            saddle_m(i, j) = 1;
            indices = [i, j; indices];
        end
    end
end
end