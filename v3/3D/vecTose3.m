function se3 = vecTose3(v)
    se3 = [skew(v(1:3)) v(4:6) ; 0 0 0 0];
end
