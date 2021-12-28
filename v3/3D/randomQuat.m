function q = randomQuat()
    q = randn(4,1);
    q = q/norm(q);
end

