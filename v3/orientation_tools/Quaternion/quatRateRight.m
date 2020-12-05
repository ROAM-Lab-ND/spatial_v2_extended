function quatDot = quatRateRight(quat,w)
    quatDot = quatProduct(quat, [0;w]/2);
end

