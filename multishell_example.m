% multishell example
[tensor,bvals] = multishell_tensor([500,1000,2000],[30,60,120],15);

output = tensor;
figure (1)
scatter3(output(:,1),output(:,2),output(:,3),25, bvals,'filled')
hold all
scatter3(-output(:,1),-output(:,2),-output(:,3),25, bvals)
axis('equal')
legend('Sampled','Unsampled')
title('Sampling Scheme')

% shells{1} = uniform_half_sphere(64,200);
% shells{2} = uniform_half_sphere(64,200);
% shells{3} = uniform_half_sphere(128,200);
% 
% nshells = optimize_shells(shells,400);
% 
% % close all
% % output = nshells{1};
% % scatter3(output(:,1),output(:,2),output(:,3),'r')
% % hold all
% % scatter3(-output(:,1),-output(:,2),-output(:,3),'r')
% % output = nshells{2};
% % scatter3(output(:,1),output(:,2),output(:,3),'g')
% % scatter3(-output(:,1),-output(:,2),-output(:,3),'g')
% % output = nshells{3};
% % scatter3(output(:,1),output(:,2),output(:,3),'b')
% % scatter3(-output(:,1),-output(:,2),-output(:,3),'b')
% % axis('equal')
% 
% tensor = [nshells{1}*sqrt(1/5); nshells{2}*sqrt(3/5); nshells{3}];
% output = tensor;
% scatter3(output(:,1),output(:,2),output(:,3))
% hold all
% scatter3(-output(:,1),-output(:,2),-output(:,3))
% axis('equal')
% 
% % randomize shell order
% p = randperm(256);
% rtensor = tensor(p,:);
% 
% % insert b-values
% i = 1;
% btensor = [];
% while i <= 256
%     btensor = [btensor; rtensor(i,:)];
%     if(mod(i,15) == 0)
%         btensor = [btensor; [0,0,0]];
%     end
%     i = i+1;
% end
% 
% bvals = sum(btensor.^2,2)*5000