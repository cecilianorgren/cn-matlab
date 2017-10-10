%% Biot-Savart integration on a generic curve
% Alessandro Masullo, 06/15/2013
clc; clear; close all
%% Domain discretization
ND = 7;
Dom = [-1.1 1;
       -1.1 1;
        0.1 6];
% Induction constant
gamma = 1;

% Integration step size
ds = 0.1;

%% Induction curve
switch menu('Choose a test case:', 'Straight wire', 'Bent wire', 'Solenoid');
    case 1
        % Test case 1: Straight wire
        L = [0 0 0;
             0 0 5];
    case 2
        % Test case 2: Bent wire
        L = [-.5 0 0;
             -.5 0 4;
             1.5 0 4];
    case 3
        % Test case 3: Solenoid
        theta = linspace(0, 15*pi, 70);
        L = [cos(theta') sin(theta') theta'/10];
end

Nl = numel(L)/3; % Number of points of the curve

%% Declaration of variables
% Induction vector components B = (U, V, W);
U = zeros(ND, ND, ND);
V = zeros(ND, ND, ND);
W = zeros(ND, ND, ND);

% Volume Mesh
[X, Y, Z] = meshgrid(linspace(Dom(1,1), Dom(1,2), ND), ...
                     linspace(Dom(2,1), Dom(2,2), ND), ...
                     linspace(Dom(3,1), Dom(3,2), ND));

%% Numerical integration of Biot-Savart law
Wait = waitbar(0, 'Integrating, please wait...');
for i = 1:ND
    for j = 1:ND
        for k = 1:ND
            waitbar(sub2ind([ND ND ND],k,j,i)/ND/ND/ND, Wait)
            % Ptest is the point of the field where we calculate induction
            pTest = [X(i,j,k) Y(i,j,k) Z(i,j,k)];
            % The curve is discretized in Nl points, we iterate on the Nl-1
            % segments. Each segment is discretized with a "ds" length step
            % to evaluate a "dB" increment of the induction "B".
            for pCurv = 1:Nl-1
                % Length of the curve element
                len = norm(L(pCurv,:) - L(pCurv+1,:));
                % Number of points for the curve-element discretization
                Npi = ceil(len/ds);
                if Npi < 3
                    close(Wait);
                    error('Integration step is too big!!')
                end
                % Curve-element discretization
                Lx = linspace(L(pCurv,1), L(pCurv+1,1), Npi);
                Ly = linspace(L(pCurv,2), L(pCurv+1,2), Npi);
                Lz = linspace(L(pCurv,3), L(pCurv+1,3), Npi);
                % Integration
                for s = 1:Npi-1
                    % Vector connecting the infinitesimal curve-element 
                    % point and field point "pTest"
                    Rx = Lx(s) - pTest(1);
                    Ry = Ly(s) - pTest(2);
                    Rz = Lz(s) - pTest(3);
                    % Infinitesimal curve-element components
                    dLx = Lx(s+1) - Lx(s);
                    dLy = Ly(s+1) - Ly(s);
                    dLz = Lz(s+1) - Lz(s);
                    % Modules
                    dL = sqrt(dLx^2 + dLy^2 + dLz^2);
                    R = sqrt(Rx^2 + Ry^2 + Rz^2);
                    % Biot-Savart
                    dU = gamma/4/pi*(dLy*Rz - dLz*Ry)/R/R/R;
                    dV = gamma/4/pi*(dLz*Rx - dLx*Rz)/R/R/R;
                    dW = gamma/4/pi*(dLx*Ry - dLy*Rx)/R/R/R;
                    % Add increment to the main field
                    U(i,j,k) = U(i,j,k) + dU;
                    V(i,j,k) = V(i,j,k) + dV;
                    W(i,j,k) = W(i,j,k) + dW;
                end
            end
        end
    end
end
close(Wait);

%% Graphic
figure(1)
M=sqrt(U.^2+V.^2+W.^2);
subplot 121
title('Normalized field')
quiver3(X,Y,Z,U./M,V./M,W./M), hold on, axis equal, grid off
plot3(L(:,1),L(:,2),L(:,3),'r-o','linewidth',3)
view(3)

subplot 122
title('Magnitude field')
quiver3(X,Y,Z,U,V,W), hold on, axis equal, grid off
plot3(L(:,1),L(:,2),L(:,3),'r-o','linewidth',3)
view(3)