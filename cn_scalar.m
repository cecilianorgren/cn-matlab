function f = cn_scalar(a,b)
sizea=size(a,2);
switch sizea
    case 3     
        f=sum(a.*b);
    case 4
        f=[a(:,1) sum(a(:,2:4).*b(:,2:4))]; 
end