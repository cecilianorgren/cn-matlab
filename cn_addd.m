function [y]=irf_add(c1,x1,c2,x2,varargin)

%IRF_ADD   add two vectors 
%
% [y]=irf_add(c1,x1,c2,x2,...,cn,xn)
% estimates y=c1*x1+c2*x2+...+cn*xn;
% c1,c2,...,cn - scalars
% x1,x2,...,xn - time series with column one being time
%
% $Id: irf_add.m,v 1.5 2010/04/14 07:59:06 andris Exp $

if mod(varargin,2)
    irf_log('proc','ERROR: Uneven number of inputs.')
    return
end

if size(varargin,2)==0
    if size(x1,2)>size(x2,2),
      irf_log('proc','WARNING: x1 has more columns than x2');
      if size(x2,2)==2, for j=3:size(x1,2), x2(:,j)=x2(:,2); end % assume x2 is time series of scalar
      elseif size(x2,2)==1, x2(:,2)=ones(size(x1))*x2(1,1);x2(:,1)=x1(:,1); % add x2(1,1) to all x1 
      elseif size(x2,2)==size(x1,2)-1, x2=[x1(1,1) x2(1,:)]; % use only first row of x2, taking time from x1
      else
        irf_log('proc','ERROR: could not make inteligent guess what you are meaning.');return; 
      end

    elseif size(x2,2)>size(x1,2),
      irf_log('proc','WARNING: x2 has more columns than x1');
      if size(x1,2)==2, for j=3:size(x2,2), x1(:,j)=x1(:,2); end % assume x2 is time series of scalar
      elseif size(x1,2)==1, x1(:,2)=ones(size(x2))*x1(1,1);x1(:,1)=x2(:,1); % add x2(1,1) to all x1 
      elseif size(x1,2)==size(x2,2)-1, x1=[x2(1,1) x1(1,:)]; % use only first row of x2, taking time from x1
      else
        irf_log('proc','irf_add() ERROR: could not make inteligent guess what you are meaning.');return; 
      end
    end

    if size(x1,1) ~= size(x2,1), 
     irf_log('proc','interpolating x2 to x1.');
     x2=irf_resamp(x2,x1);
    end

    y=x1;

    y(:,2:end)=c1*x1(:,2:end)+c2*x2(:,2:end);
    
elseif size(varargin,2)==2
    y=cn_addd(1,cn_addd(c1,x1,c2,x2),varargin{1},varargin{2});
else
    y=cn_addd(1,cn_addd(c1,x1,c2,x2),varargin{1},varargin{2},varargin{3:end});
end
