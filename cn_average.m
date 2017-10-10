function f=cn_average(a,varargin)
args=size(varargin,1);
switch args
    case 0 % sum over rows
        a_cols=size(a,2);
        a_rows=size(a,1);
        switch a_cols
            case 2 % for example temperature or density
                av=sum(a(:,2)/a_rows,1);
                f=av;
            case 3
                av=sum(a,1)/a_rows;
                f=av;
            case 4
                if a_rows==1
                    index=1;
                else
                    index=fix(a_rows/2);
                end
                av=sum(a(:,2:4),1)/a_rows;       
                f=[a(index,1) av];
            case 5
                if a_rows==1
                    index=1;
                else
                    index=fix(a_rows/2);
                end
                av=sum(a(:,2:5),1)/a_rows;      
                f=[a(index,1) av];           
        end 
    case 1 % only two row vectors
        b=varargin{1,1};
        a_cols=size(a,2);
        a_rows=size(a,1);
        switch a_cols
            case 2
                av=(a(:,2)+b(:,2))/2;  
                
                if a_rows==1
                    index=1;
                else
                    index=fix(a_rows/2);
                end
                     
                f=[a(1) av]; 
            case 3
                av=(a+b)/a;
                f=av;
            case 4
                av=(a(2:4)+b(2:4))/2;       
                f=[a(1) av];
            case 5
                av=(a(2:5)+b(2:5))/2;     
                f=[a(1) av];           
        end 
end