function tens_out = tensorRotT(tens_in)
% takes one by one transpose of tensor - along 1-2 dims

tens_out = permute(tens_in,[2 1 3]);