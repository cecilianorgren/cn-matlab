function B = cn_m_trans(A,M,way)
% cn_m_trans(time vector,M,flag)
%
%flag=1 : GSE to field aligned
%flag=-1 : field aligned to GSE

if size(A,2)>3
    si=2;
    B=A(:,1:4);
else
    si=1;
    B=A;
end


switch way
    case 1  % from gse to m    
        for k=1:size(A,1)
            B(k,si:si+2)=(M*A(k,si:si+2)')';
        end
    case -1 % from m to gse
        for k=1:size(A,1)
            B(k,si:si+2)=(M\A(k,si:si+2)')';
        end
end