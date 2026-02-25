function f = cn_mag(a)
sizea=size(a,2);
switch sizea
    case 1 % column vector
        f=sqrt( a(1).^2 + a(2).^2 + a(3).^2 );
    case 3
        f=sqrt( a(:,1).^2 + a(:,2).^2 + a(:,3).^2 );
    case 4
        f=[a(:,1) sqrt( a(:,2).^2 + a(:,3).^2 + a(:,4).^2 )];
    case 5 % 5th col of B is already abs
        f=[a(:,1) a(:,5)];   
        f=[a(:,1) sqrt( a(:,2).^2 + a(:,3).^2 + a(:,4).^2 )];
end    