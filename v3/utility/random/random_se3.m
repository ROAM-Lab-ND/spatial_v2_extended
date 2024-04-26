function T = random_se3()
    p = randn(3,1);
    W = randomSkew(3);
    T = [W p ; 0 0 0 0];
end