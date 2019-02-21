function split_directions(tensor,i)
% prints a tensor file with diffusion scheme spec. by Nx3 vector tensor across multiple acquisitions.  i is a vector
% containing the number of directions in each acquisition.
%i = 77:84;
counter = 1;
for ii = 1:length(i)
    fprintf([num2str(i(ii)) '\n'])
    for jj = 1:i(ii)
        fprintf([num2str(tensor(counter,:)) '\n'])
        counter = counter + 1;
    end
        
end