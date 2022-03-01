function tens_out = rotR(tens_in)

m = size(tens_in,3);
n = size(tens_in,2);

 for jj=1:m
   for ii=1:n
        temp1 = tens_in(:,ii,jj);
        tens_out(:,jj,ii) = temp1;
    end    
end