function j = cn_j(deltab,distance,option)
% CN_J(deltaB,distance)
%
% distance in km (is resampled to timline of delta B)
% delta B in nT
% returns current in A

[row col]=size(distance);

Units=irf_units;

if col==3;
    distance=irf_resamp([1 distance],deltab);
else
    distance=irf_resamp(distance,deltab,'linear');
end

deltab(:,2:4)=deltab(:,2:4)*1e-9;       % b in T
distance(:,2:4)=distance(:,2:4)*1e3;    % distance in m
nabla=distance;
nabla(:,2:4)=nabla(:,2:4).^(-1);        % 1/distance

j=irf_cross(nabla,deltab);              % nabla x b = mu0*j

j(:,2:4)=j(:,2:4)/Units.mu0; % T m gives A

