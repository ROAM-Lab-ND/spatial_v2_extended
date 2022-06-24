function Y_group = groupRegressor(a,v)
    N = length(v)/6;
    
    a = reshape(a, [6 N]);
    v = reshape(v, [6 N]);
    
    group_inds = @(i,n) ((i-1)*n+1):(n*i);
    rows = @(i) group_inds(i,6);
    cols = @(i) group_inds(i,10);
    
    Y_group = zeros(6*N, 10*N);
    for i = 1:N
        Y_group( rows(i), cols(i) ) = individualRegressor(a(:,i), v(:,i));
    end
end