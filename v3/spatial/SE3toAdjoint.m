function X = SE3toAdjoint(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    O = zeros(3);
    X = [R O; skew2(p)*R R];
end
