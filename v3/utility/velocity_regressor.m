function y_V = velocity_regressor(V)
% assemble velocity regressor based on twist V

y_V = zeros(6,10);
w = V(1:3); v = V(4:6);
y_V(4:6,1) = v;
y_V(:,2:4) = [-skew(v);skew(w)];
y_V(1:3,5:end) = [w(1),  0  ,  0   ,  w(2),   0  , w(3);
                               0  , w(2),  0   ,  w(1),  w(3),   0 ;
                               0  ,  0  , w(3) ,   0  ,  w(2), w(1) ];

end