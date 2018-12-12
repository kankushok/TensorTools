function[output] = uniform_half_sphere(num_points,iter)

%this program generates points evenly distributed over the
%unit sphere by modeling electrostatic repulsion.  The inputs are the
%number of
tic
%num_points = 400;

%plot progress
display = 0;
hist = zeros(1,iter);
%generate a rotation matrix which rotates
rotMtx = @(u,theta)[cos(theta)+u(1).^2*(1-cos(theta)), u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta), u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta);...
    u(2)*u(1)*(1-cos(theta))+u(3)*sin(theta), cos(theta)+u(2).^2*(1-cos(theta)), u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta); ...
    u(3)*u(1)*(1-cos(theta))-u(2)*sin(theta), u(3)*u(2)*(1-cos(theta))+u(1)*sin(theta) cos(theta)+u(3).^2*(1-cos(theta))];

% generate num_points on a sphere randomly
allpts = randn(num_points,3);
allpts = allpts./repmat(sqrt(sum(allpts.^2,2)),1,3);
allpts = [allpts; -allpts];
for ii = 1:iter
    % generate force vectors on each point
    F = zeros(num_points*2,num_points*2,3);
    F(:,:,1) = bsxfun(@minus,allpts(:,1),allpts(:,1)');
    F(:,:,2) = bsxfun(@minus,allpts(:,2),allpts(:,2)');
    F(:,:,3) = bsxfun(@minus,allpts(:,3),allpts(:,3)');
    F = F./repmat(sum(F.^2,3).^1.5,1,1,3);
    F(isnan(F)) = 0;
    netF = squeeze(sum(F,2));

    % reject Force along i
    netF = netF - repmat(sum(netF.*allpts,2),1,3).*allpts;
    
    % calculate torque
    torque = bsxfun(@cross,allpts,netF);
    
    % calculate net torque
    torque = repmat(torque(1:end/2,:) + torque(end/2+1:end,:),2,1,1);
    % calculate magnitude (rotation angle)
    mtorque = sqrt(sum(torque.^2,2));
    
    % normalized rotation axes
    ntorque = bsxfun(@rdivide, torque, mtorque);
    
    % set max step size to 1 degree
    scale = pi/180/max(abs(mtorque));
    if scale > 1e6
        return;
    end
    
    % update positions
    for jj = 1:num_points*2
        allpts(jj,:) = (rotMtx(ntorque(jj,:),mtorque(jj)*scale)*allpts(jj,:)')';
    end
    
    % debugtool
    if(display == 1)
        hold off
        scatter3(allpts(:,1),allpts(:,2),allpts(:,3),'r')
        axis('equal')
        hold on
        quiver3(allpts(:,1),allpts(:,2),allpts(:,3),netF(:,1)/10,netF(:,2)/10,netF(:,3)/10,0);
        axis([-1 1 -1 1 -1 1])
        pause
    end
end

%output = allpts(1:num_points,:);

% pick half directions
output = allpts(1,:);
remainder = allpts([2:end/2,end/2+2:end],:);
%draft by largest distance to (previously) chosen points
for ii = 1:num_points-1
    clear d
    d(:,:,1) = bsxfun(@minus,remainder(:,1),output(:,1)');
    d(:,:,2) = bsxfun(@minus,remainder(:,2),output(:,2)');
    d(:,:,3) = bsxfun(@minus,remainder(:,3),output(:,3)');
    d = sum(d.^2,3); % squared distance
    d = mean(d,2); % mean squared distance compared to output points
    idx = find(d==max(d));
    N = num_points - ii;
    keep = ones(1,2*N);
    keep(idx) = 0;
    if(idx <= N)
        keep(idx+N) = 0;
    else
        keep(idx-N) = 0;
    end
    % update remainder and output
    output = [output; remainder(idx,:)];
    remainder = remainder(keep==1,:);
end


toc