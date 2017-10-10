function f=cn_toepoch(tint,varargin)

if size(tint)==1
    ep_tint=tint;
else
    ep_tint=toepoch(tint(1:6));%+tint(7)*0.01;
end

switch size(varargin,2)
    
    case 0 % If only one parameter is given 
        
        f=ep_tint;
        
    case 1 % A time series is also given 
        
        A=varargin{1,1};
        A_index1=find(A(:,1)<ep_tint,1,'last');
        A_index2=find(A(:,1)>ep_tint,1,'first');
        
        A_t1=A(A_index1,:); % closest over
        A_t2=A(A_index2,:); % closest under
        
        % Check which one is closer in time
        
        d1=A_t1(1)-ep_tint;
        d2=ep_tint-A_t2(1);
        
        if (d1>d2)
            f=A_t2;
        elseif (d2>d1)
            f=A_t1; 
        else
            f=A_t1; 
        end
        
    case 2 % Two time + a time series 
        
        A=varargin{1,2};
        tint2=varargin{1,1};
        if size(tint)==1
            ep_tint2=tint2;
        else
            ep_tint2=toepoch(tint2(1:6));%+tint(7)*0.01;
        end
        %ep_tint2=toepoch(tint2);%(1:6))+tint2(7)*0.01;
        A_index1=find(A(:,1)<ep_tint,1,'last');
        A_index2=find(A(:,1)>ep_tint2,1,'first');
        
        A_t12=A(A_index1:A_index2,:);
        f=A_t12;
        
end 