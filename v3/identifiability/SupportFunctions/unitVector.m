function e = unitVector(i,n)
%% e = unitVector(i,n)
% gives the unit vector such that length(e) = n and e_i = 1
    e = zeros(n,1);
    e(i) = 1;
end
