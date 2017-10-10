function [av hat mag]=cn_hat(a,varargin)
args=size(varargin,1);
a_cols=size(a,2);
a_rows=size(a,1);

switch args
    case 0 % only one vector
        aav=sum(a(:,2:4),1)/a_rows;
        av=aav;
        mag=sqrt(aav(1)^2+aav(2)^2+aav(3)^2);
        hat=aav/mag;
    case 1 % two vectors
        b=varargin{1,1};
        b_cols=size(b,2);
        b_rows=size(b,1);
        
        aav=sum(a(:,2:4),1)/a_rows;
        bav=sum(b(:,2:4),1)/b_rows;
        av=(aav+bav)/2;
        mag=sqrt(av(1)^2+av(2)^2+av(3)^2);
        hat=av/mag;
end
end