function T = AdjointToSE3(X)
    R = X(1:3,1:3);
    p = skew( X(4:6,1:3)*R' );
    T = [R p; 0 0 0 1];
end
