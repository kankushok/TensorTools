function [ btensor, bvals ] = multishell_tensor( bvals, numdirs, num_b0 )
%multishell_tensor This function generates the gradient directions suitable for use with the tensor.dat file
%read by epi2. It automatically generates uniformly sampled
%directions on a half sphere for each shell and optimizes the rotation
%between shells for optimal angular resolution
%   bvals is 1xNshells vector containing the b-value for each shell
%   numdirs is 1xNshells containing the number of directions on each shell
%   num_b0 < sum(numdirs) is the number of b0s interspersed through the diffusion volume.

Nshells = length(bvals);
Ndiff = sum(numdirs);

for ii = 1:Nshells
    shells{ii} = uniform_half_sphere(numdirs(ii),200);
end
if(Nshells > 1)
    nshells = optimize_shells(shells,400);
else
    nshells = shells;
end
tensor = [];
for ii = 1:Nshells
    tensor = [tensor; nshells{ii}*sqrt(bvals(ii)/max(bvals))]; 
end

% randomize shell order
p = randperm(Ndiff);
rtensor = tensor(p,:);

% insert b-values
i = 1;
btensor = [];
interval = floor(Ndiff/num_b0);
while i <= length(tensor)
    btensor = [btensor; rtensor(i,:)];
    if(mod(i,interval) == 0)
        btensor = [btensor; [0,0,0]];
    end
    i = i+1;
end

bvals = sum(btensor.^2,2)*max(bvals);
end

