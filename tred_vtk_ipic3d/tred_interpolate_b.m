function results = tred_interpolate_b(core,satellite,timestep,thickness)
%
%
%   Finds to what core the B-field points.
%



% Load B
eval(['fid=fopen(''/Users/Cecilia/TRED46/VirtualSatellite/VirtualSatelliteTraces',num2str(core),'.txt'');'])
trace_format = repmat('%f',1,14*27);
data=textscan(fid, trace_format, 1, 'headerlines', 27+timestep-1);
b=[data{1,satellite*14+1} data{1,satellite*14+2} data{1,satellite*14+3}];

% Find all positions
csp=tred_find_cores([0 20],[0 15],[0 10]);
[~,index]=ismember([core satellite],[csp{1} csp{2}],'rows');

% Take some close to B-direction
surrounding_satellites=tred_surr_sat(core,satellite,thickness);

% Find distance from satellites to B-line
p1=csp{3}(index,:);
p2=p1+irf_norm(b);
p2_minus=p1-irf_norm(b);

for k=1:((thickness*2+1)^3)
    [~,index]=ismember([surrounding_satellites{1}(k) surrounding_satellites{2}(k)],[csp{1} csp{2}],'rows');
    p3(k,:)=csp{3}(index,:);
    u(k)=(p3(k,1)-p1(1))*(p2(1)-p1(1))+(p3(k,2)-p1(2))*(p2(2)-p1(2))+(p3(k,3)-p1(3))*(p2(3)-p1(3));
    x(k)=p1(1)+u(k)*(p2(1)-p1(1));
    y(k)=p1(2)+u(k)*(p2(2)-p1(2));
    z(k)=p1(3)+u(k)*(p2(3)-p1(3));
    distance(k,1)=sqrt((x(k)-p3(k,1))^2+(y(k)-p3(k,2))^2+(z(k)-p3(k,3))^2);
end

1;
plot3([p2_minus(1) p2(1)],[p2_minus(2) p2(2)],[p2_minus(3) p2(3)],'r'); hold on
plot3(p1(1),p1(2),p1(3),'go')
for k=1:((thickness*2+1)^3)
    if distance(k,1)<0.5; color='g'; else; color='b'; end;
    plot3([p3(k,1) x(k)],[p3(k,2) y(k)],[p3(k,3) z(k)],color); hold on
end

axis equal

results={surrounding_satellites{1},surrounding_satellites{2},...
         surrounding_satellites{3},distance};