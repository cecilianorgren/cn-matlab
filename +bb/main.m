% program to show bouncing balls/bubbles.
% set colormap
whitered = [];
for ii=1:100
    whitered = [whitered; [1 1*ii/100 1*ii/100]];
end
colormap(whitered)
% array with structures containing all bubble properties
maxr = 1;
maxv = 0.1;
maxn = 1;
x = [0 1]*10;
y = [0 1]*10; 

nBubbles = 30;
vecBubbles = cell(nBubbles,1);
totBubbles = 0;
while totBubbles<nBubbles
    totBubbles = totBubbles + 1;
    vAngle = 2*pi*rand(1,1);
    bubble.vx = maxv*cos(vAngle);
    bubble.vy = maxv*sin(vAngle);
    bubble.r = rand(1,1)*maxr;
    bubble.x = x(1)+bubble.r + (diff(x)-2*bubble.r)*rand(1,1);
    bubble.y = y(1)+bubble.r + (diff(y)-2*bubble.r)*rand(1,1);
    bubble.n = maxn*rand(1,1);
    bubble.m = bubble.r^2*pi*bubble.n;
    vecBubbles{totBubbles} = bubble;    
    
    % check if they overlap
    for ii = 1:(totBubbles-1)
        realDistance = sqrt((vecBubbles{ii}.x-bubble.x)^2+(vecBubbles{ii}.y-bubble.y)^2);
        minDistance = vecBubbles{ii}.r + bubble.r;
        if realDistance<minDistance
            totBubbles = totBubbles - 1;
            %disp('Bubbles overlap. Trying again.')
        end   
    end
end

bb.plot(vecBubbles,x,y);


% time stepping
nt = 100;
tstop = 5;
dt = tstop/nt;
t = 0;

while t<tstop;
    % circle through all the particles to see if they will touch and update
    % position
    % if they will touch each other or the wall, change the velocity vector
    
    for ii = 1:(nBubbles-1)
        % update position
        vecBubbles{ii}.x = vecBubbles{ii}.x + dt*vecBubbles{ii}.vx;
        vecBubbles{ii}.y = vecBubbles{ii}.y + dt*vecBubbles{ii}.vy;
        
        % check if bubble now touches wall
        if vecBubbles{ii}.x + vecBubbles{ii}.r > x(2) && vecBubbles{ii}.vx > 0;
            vecBubbles{ii}.vx = -vecBubbles{ii}.vx;
        elseif vecBubbles{ii}.x - vecBubbles{ii}.r < x(1) && vecBubbles{ii}.vx < 0;
            vecBubbles{ii}.vx = -vecBubbles{ii}.vx;
        elseif vecBubbles{ii}.y + vecBubbles{ii}.r > y(2) && vecBubbles{ii}.vy > 0;
            vecBubbles{ii}.vy = -vecBubbles{ii}.vy;
        elseif vecBubbles{ii}.y - vecBubbles{ii}.r < y(1) && vecBubbles{ii}.vy < 0;
            vecBubbles{ii}.vy = -vecBubbles{ii}.vy;
        end
        % check if bubble collide with other bubble
        for jj = (ii+1):nBubbles
            %disp(num2str([ii jj])) 
            % vector from ii to jj
            sepVector = [vecBubbles{jj}.x-vecBubbles{ii}.x vecBubbles{jj}.y-vecBubbles{ii}.y];
            realDistance = sqrt(sepVector(1).^2+sepVector(2).^2);
            %realDistance = sqrt((vecBubbles{ii}.x-vecBubbles{jj}.x)^2+(vecBubbles{ii}.y-vecBubbles{jj}.y)^2);
            minDistance = vecBubbles{ii}.r + vecBubbles{jj}.r;
            vdot = dot([vecBubbles{ii}.vx vecBubbles{ii}.vy],[vecBubbles{jj}.vx vecBubbles{jj}.vy]);            
            if realDistance<minDistance
                %disp(vdot)
                if vdot > 0; continue; end
                %disp('Bubbles touch. Redirecting.')
                % first particle
                v = [vecBubbles{ii}.vx vecBubbles{ii}.vy];
                normv = dot(sepVector,v);
                newv = v - 2*sepVector/norm(sepVector)*normv;
                vecBubbles{ii}.vx = newv(1);
                vecBubbles{ii}.vy = newv(2);
                %disp(num2str([v norm(v)]))
                %disp(num2str([newv norm(newv)]))
                %tempA = sepVector;
                %tempBii = [vecBubbles{ii}.vx vecBubbles{ii}.vy];
                %inpactAngleii = acos(dot(tempA/norm(tempA),tempBii/norm(tempBii)));
                
                % second particle                
                v = [vecBubbles{jj}.vx vecBubbles{jj}.vy];
                normv = dot(sepVector,v);
                newv = v - 2*sepVector/norm(sepVector)*normv;
                vecBubbles{jj}.vx = newv(1);
                vecBubbles{jj}.vy = newv(2);
                %tempA = sepVector;
                %tempBjj = [vecBubbles{jj}.vx vecBubbles{jj}.vy];
                %inpactAnglejj = acos(dot(tempA/norm(tempA),tempBjj/norm(tempBjj)));
                
                %disp(num2str([inpactAngleii inpactAnglejj tempBii tempBjj]))                
            end   
            
        end
    end
    % plot progress
    bb.plot(vecBubbles,x,y);
    pause(0.1)
    t = t + dt;
end