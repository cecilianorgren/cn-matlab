function tred_plot_pos(cores,sc,Position)
%
%
%
axis([floor(min(Position(:,1))) ceil(max(Position(:,1))),...
      floor(min(Position(:,2))) ceil(max(Position(:,2))),...
      floor(min(Position(:,3))) ceil(max(Position(:,2)))])
npos=size(cores,1);

for k=1:npos       
    cor=cores(k); % cores
    sat=sc(k); % satellites      
    %quiver3(Position(k,sat,1),Position(k,sat,2),Position(k,sat,3),...
    %        Bnorm(k,sat,16000,1)*0.3,Bnorm(k,sat,16000,2)*0.3,Bnorm(k,sat,16000,3)*0.3); hold on;        
    text(Position(k,1),Position(k,2),Position(k,3),[num2str(cor),'-',num2str(sat)]);
end
xlabel('x: xgsm');ylabel('y: zgms');zlabel('z: -ygsm')