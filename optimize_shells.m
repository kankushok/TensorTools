function[nshells] = optimize_shells(shells,iter)

%this program generates rotates b-shells to maximize angular resolution by
%modeling by modeling electrostatic repulsion.  The inputs are a cell
%containing the Nx3 matrix describing the directions in that shell. Iter is
%number of iterations to run
%
tic

num_shells = length(shells); 
allpts = [];
shellid = [];
for i = 1:num_shells
    num_pts(i) = size(shells{i},1);
    allpts = [allpts ; shells{i}];
    shellid = [shellid ; i*ones(num_pts(i),1)];
end

allpts = [allpts; -allpts];
shellid = [shellid; shellid];
N = size(allpts,1);
%plot progress
display = 0;

%generate a rotation matrix which rotates
rotMtx = @(u,theta)[cos(theta)+u(1).^2*(1-cos(theta)), u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta), u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta);...
    u(2)*u(1)*(1-cos(theta))+u(3)*sin(theta), cos(theta)+u(2).^2*(1-cos(theta)), u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta); ...
    u(3)*u(1)*(1-cos(theta))-u(2)*sin(theta), u(3)*u(2)*(1-cos(theta))+u(1)*sin(theta) cos(theta)+u(3).^2*(1-cos(theta))];

for ii = 1:iter
    % generate force vectors on each point
    F = zeros(N,N,3);
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
    for jj = 1:num_shells
        ntorque = sum(torque(shellid==jj,:),1);
        % calculate magnitude (rotation angle)
        mtorque = sqrt(sum(ntorque.^2,2));
        
        % normalized rotation axes
        %ntorque = bsxfun(@rdivide, ntorque, mtorque);
        ntorque = ntorque/mtorque;
        
        % set max step size to 1 degree
        scale = pi/180/max(abs(mtorque));
        if scale > 1e10
            return;
        end
        
        % update positions
        allpts(shellid==jj,:) = (rotMtx(ntorque,mtorque*scale)*allpts(shellid==jj,:)')';
        
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
end

for i = 1:num_shells
    tmp = allpts(shellid==i,:);
    nshells{i} = tmp(1:end/2,:);
end


toc