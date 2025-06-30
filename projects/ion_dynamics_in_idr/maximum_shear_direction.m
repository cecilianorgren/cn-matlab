function varargout = maximum_shear_direction(P)

T = P.data;

% Rotate tensor, to find the most unequal component
theta1 = 0.5*atan((2*T(:,2,3))./(T(:,3,3)-T(:,2,2)));
theta2 = 0.5*atan((2*T(:,1,2))./(T(:,2,2)-T(:,1,1)));
%theta2 = 0.5*atan((2*T(:,2,1))./(T(:,1,1)-T(:,2,2)));

% Initialize arrays
R2 = zeros(P.length,3,3);
T2 = zeros(P.length,3,3);

R2(:,1,:) = [cos(theta2) -sin(theta2) theta2*0];
R2(:,2,:) = [sin(theta2) cos(theta2) theta2*0];
R2(:,3,:) = repmat([0 0 1],P.length,1);

%T2(:,1,:) = R2(:,1,:)*(T(:,:,1).*R2(:,1) + T(:,:,2).*R2(:,2) + T(:,:,3).*R2(:,3));


for it = 1:P.length
  R2tmp = squeeze(R2(it,:,:));
  Ttmp = squeeze(T(it,:,:));
  T2tmp = R2tmp*(Ttmp)*transpose(R2tmp);
  
  if T2tmp(1,1) > T2tmp(2,2)
    theta2(it,:) = -theta2(it,:,:);
    R2tmp = [cos(theta2(it,:)) -sin(theta2(it,:)) 0; sin(theta2(it,:)) cos(theta2(it,:)) 0; 0 0 1];
    T2tmp = R2tmp*(Ttmp*transpose(R2tmp));
  end
  
  T2(it,:,:) = T2tmp;
end


switch nargout
  case 1
    varargout{1} = T2;
  case 2
    varargout{1} = T2;
    varargout{2} = R2;
end
