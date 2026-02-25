function plot(bubbles,x,y)
% bubbles is cell array of bubble structures

%h = gca;
%hold(h,'off')
oldObj = findobj(gca,'type','patch');
delete(oldObj);
nBubbles = numel(bubbles);
angle = linspace(0,2*pi,90);
axis equal
for ii = 1:nBubbles
    cx = bubbles{ii}.x + bubbles{ii}.r*cos(angle);
    cy = bubbles{ii}.y + bubbles{ii}.r*sin(angle);
    patch(cx,cy,bubbles{ii}.n);
    set(gca,'xlim',x,'ylim',y);
    hold(gca,'on');
end
hold(gca,'off');