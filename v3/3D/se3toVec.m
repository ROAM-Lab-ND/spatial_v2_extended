function v = se3toVec(se3)
    v = [skew(se3(1:3,1:3)) ; se3(1:3,4)];
end
