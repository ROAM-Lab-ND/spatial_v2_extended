function X = randomAdSE3()
    p = randn(3,1);
    R = randomRotation();
    X = plux(R,p);
end
