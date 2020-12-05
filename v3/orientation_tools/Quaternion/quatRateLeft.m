function quatDot = quatRateLeft(quat,w)
    quatDot = quatProduct([0;w]/2, quat);
end

