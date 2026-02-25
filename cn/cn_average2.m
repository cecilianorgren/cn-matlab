function aav=cn_average2(a,how)
a_cols=size(a,2);
a_rows=size(a,1);
switch how
    case 1 % sum over rows
        aav=sum(a,1)/a_rows;
        aav=[ aav(:,2:a_cols)];
    case 2 % sum over columns
        aav=[a(:,1) sum(a(:,2:a_cols),2)];
    case 3
        aav1=sum(a,1)/a_rows;
        aav=[a(:,1) sum(aav1(:,2:a_cols),2)];
end
end