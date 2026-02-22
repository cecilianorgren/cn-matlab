function polys = splitContours(C)
    k = 1;
    polys = {};
    while k < size(C,2)
        level = C(1,k);
        npts  = C(2,k);
        idx = k+1 : k+npts;

        x = C(1,idx);
        y = C(2,idx);

        polys{end+1} = struct('level', level, 'x', x, 'y', y);

        k = k + npts + 1;
    end
end
