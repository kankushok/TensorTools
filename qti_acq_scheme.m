function [waveform,dir] = qti_acq_scheme(sde_bvals, ide_bvals, num_b0, num_sde, num_ide)
%qti_acq_scheme generates tensor and waveform files for qti acq. with 
%   inputs:
%   sde_bvals vector of bvals for sde acq.
%   ide_bvals vector of bvals for ide acq.
%   num_b0 scalar number of b0s to interleave into acq.
%   num_sde vector with number of acq. per bval sde
%   num_ide vector with number of acq. per bval ide

% sde_bvals = [ 1000,1500,2000]
% num_sde = [30,30,60]
% ide_bvals = [1000,1250,1500,2000]
% num_ide = [8,8,8,56]
% num_b0 = 10;

% set up SDE acquisitions
[tensor,bvals] = multishell_tensor(sde_bvals,num_sde,0);


% isotropic encoding 
qte = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
qte = repmat(qte,ceil(sum(num_ide)/8), 1);

scale_ide = [];
for i = 1:length(ide_bvals)
    scale_ide = [scale_ide ; ones(num_ide(i),3)*sqrt(ide_bvals(i)/max(ide_bvals))];
end

dir = [tensor; qte.*scale_ide];

% tells sequence which waveforms to play out
waveform = ones(sum(num_sde) + sum(num_ide),1);
waveform(1:sum(num_sde)) = 1;
waveform(sum(num_sde)+1:end) = 3;

% scrambles waveform order randomly...
k = randperm(size(dir,1));
dir = dir(k,:);
waveform = waveform(k);

% insert b0s at regular intervals...
k = floor(size(dir,1)/num_b0);
for i = k:k:length(dir)
    dir(i+2:end+1,:) = dir(i+1:end,:);
    waveform(i+2:end+1) = waveform(i+1:end);
    dir(i+1,:) = [0 0 0];
    waveform(i+1) = 0;
end

end
