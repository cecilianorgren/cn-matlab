function out = smooth2d(A,np)

[n1,n2] = size(A);
B = zeros(n1,n2);
ii2 = 1:n2;

% Wrong
for i2 = 1:n2  
  % First smooth over first dimension
  B(:,i2) = smooth(A(:,i2),2*np+1,'moving'); % overwrite boundary values
  Ashift = circshift(A,2*np-1,1);
  Btmp = Ashift;
  Btmp_smooth = smooth(Ashift(:,i2),2*np+1,'moving');
  Btmp_smooth_shiftback = circshift(Btmp_smooth,-(2*np-1),1);
  
  B(1:np,i2) = Btmp_smooth_shiftback(1:np);
  B(n1-np+1:n1,i2) = Btmp_smooth_shiftback(n1-np+1:n1);
  %B(n1-np:n1,i1) = Btmp(np+np+1:np+2*np,i1);
  %plot(ii1,A(:,i1),ii1,Ashift(:,i1),ii1,Btmp(:,i1),'--',ii1,B(:,i1),'--')
  %plot(ii2,A(:,i2),ii2,Ashift(:,i2),ii2,Btmp_smooth,'--',ii2,Btmp_smooth_shiftback,'--',ii2,B(:,i2),'--')
  %legend({'A','Ashift','Btmp smooth','Btmp smooth shift back','B'},'Location','best')  
end

for i1 = 1:n1
  % Then smooth over second dimension
  B(i1,:) = smooth(A(i1,:),2*np+1,'moving'); % overwrite boundary values
  Ashift = circshift(A,2*np-1,2);
  Btmp = Ashift;
  Btmp_smooth = smooth(Ashift(i1,:),2*np+1,'moving');
  Btmp_smooth_shiftback = circshift(Btmp_smooth,-(2*np-1),2);
  
  B(i1,1:np) = Btmp_smooth_shiftback(1:np);
  B(i1,n2-np+1:n2) = Btmp_smooth_shiftback(n2-np+1:n2);
end

out = B;