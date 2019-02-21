function[allpts] = uniform_sphere(num_points,iter)

%this program generates points evenly distributed over the
%unit sphere by modeling electrostatic repulsion.  The inputs are the
%number of
tic

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
for ii = 1:iter
    % generate force vectors on each point
    F = zeros(num_points,num_points,3);
    F(:,:,1) = bsxfun(@minus,allpts(:,1),allpts(:,1)');
    F(:,:,2) = bsxfun(@minus,allpts(:,2),allpts(:,2)');
    F(:,:,3) = bsxfun(@minus,allpts(:,3),allpts(:,3)');
    F = F./repmat(sum(F.^2,3).^1.5,1,1,3);
%     for i = 1:num_points
%         for j = 1:num_points
%             % force on point i from point j
%             F(i,j,:) = 1./norm(allpts(i,:)-allpts(j,:)).^3.*(allpts(i,:)-allpts(j,:));
%         end
%     end
    F(isnan(F)) = 0;
    netF = squeeze(sum(F,2));
    
    % reject Force along i
    netF = netF - repmat(sum(netF.*allpts,2),1,3).*allpts;
    
    magF = sum(netF,2);
    thresh = 1000;
    nmagF = magF;
    nmagF(magF>thresh) = 1000;
    
    % renormalize netF
    netF = netF./repmat(magF./nmagF,1,3);
    
    % update positions
    allpts = allpts + 0.1*netF/num_points;
    
    % renormalize allpts
    allpts = allpts./repmat(sqrt(sum(allpts.^2,2)),1,3);
    if(display==1)
        hold off
        scatter3(allpts(:,1),allpts(:,2),allpts(:,3),'r')
        axis('equal')
        hold on
        quiver3(allpts(:,1),allpts(:,2),allpts(:,3),netF(:,1),netF(:,2),netF(:,3),0);
        axis([-1 1 -1 1 -1 1])
        pause
    end
end
toc