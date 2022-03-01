function tens_out = rotT(tens_in)
% takes one by one transpose of tensor - along 1-2 dims

for ii=1:size(tens_in,3)
    tens_out(:,:,ii) = tens_in(:,:,ii).';
end