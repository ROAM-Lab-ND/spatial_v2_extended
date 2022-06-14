function X = SE3toAdjoint(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    O = zeros(3);
    X = [R O; skew(p)*R R];
end
