function T = randomSE3()
    p = randn(3,1);
    R = randomRotation();
    T = [R p ; 0 0 0 1];
end
